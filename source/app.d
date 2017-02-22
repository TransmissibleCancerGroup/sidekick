import bio.bam.read: BamRead;
import bio.bam.region: BamRegion;
import bio.bam.reader: BamReader;
import bio.bam.multireader: MultiBamReader;
import bio.bam.writer: BamWriter;
import bio.bam.pileup: pileupChunks;
import core.memory;
import std.algorithm: filter, reduce, sort;
import std.algorithm.iteration: map;
import std.container.rbtree: redBlackTree;
import std.conv: to;
import std.file: exists, remove, rmdirRecurse, mkdirRecurse;
import std.getopt;
import std.parallelism: TaskPool;
import std.path: buildPath;
import std.range: take, tee;
import std.stdio: writeln, writefln;
import utils.progressbar: ProgressBar;

/**
  Returns true for reads in read pairs with at least one unmapped read
*/
bool passes_initial_checks(R)(const R r) {
    return (!r.is_duplicate &&
            r.is_paired &&
            !r.proper_pair &&
            !r.failed_quality_control &&
            !r.is_supplementary &&
            !r.is_secondary_alignment &&
            (r.is_unmapped || r.mate_is_unmapped));
}

/**
  Returns true for a mapped read that passes a mapping quality threshold
  or an unmapped read that passes a per-base quality threshold.
*/
bool passes_quality_checks(R)(const R r, int base_qual, int map_qual) {
    return r.is_unmapped ? avg_base_quality(r) >= base_qual : r.mapping_quality >= map_qual;
}

/**
  Compute the average per-base read quality
*/
double avg_base_quality(R)(const ref R r) {
    auto qual = reduce!"a + b"(0.0, r.base_qualities);
    return qual / r.base_qualities.length;
}

/**
  Returns true if a read and a region overlap
*/
bool overlaps(Read, Region)(Read read, Region region, bool use_mate_info=false) {
    auto read_start = use_mate_info ? read.mate_position : read.position;
    auto read_ref_id = use_mate_info ? read.mate_ref_id : read.ref_id;
    return (read_ref_id == region.ref_id &&
            read_start <= region.end &&
            read_start + (use_mate_info ? read.sequence_length : read.basesCovered()) > region.start);
}

/**
  Opens a new Bam file for writing, with header info copied from
  a template reader

  The writer will need to have its .finish method called.
*/
auto new_bamwriter_from_template(Reader)(string filename, ref Reader reader, ref TaskPool taskpool) {
    // Remember to call scope(exit) writer.finish() on the returned writer
    auto writer = new BamWriter(filename, -1, taskpool);
    writer.writeSamHeader(reader.header);
    writer.writeReferenceSequenceInfo(reader.reference_sequences);
    return writer;
}

/**
  Search subject file for pairs present in query file,
  and write matching reads to output file. Uses tmpdir
  to write intermediate results before merging.
*/
void filter_bam(string query, string subject, ref TaskPool taskpool, string tmpdir, string outfile, int at_a_time=1000000) {
    import core.memory;

    string[] tmpfiles;
    auto query_reader = new BamReader(query, taskpool);

    auto queries = query_reader.reads();
    int batch = 0;
    auto cache = redBlackTree!string(); // could use TTree or HashSet from emsi_containers
    while(!queries.empty) {             // but seems there is no memory advantage
        
        scope(exit) {
            batch++;
            cache.clear();
            GC.collect();
            GC.minimize();
        }

        // populate the (empty) cache
        assert(cache.empty);
        cache.insert(queries.take(at_a_time).map!(r => r.name));
        if (cache.length == 1) break;  // maybe queries.empty is broken? 
        writefln("[filter_bam] - batch number %d processing %d reads", batch, cache.length);

        // Open subject file
        auto subject_reader = new BamReader(subject, taskpool);

        // Open writer for this iteration
        string tmpfilename = buildPath(tmpdir, "tmp" ~ to!string(batch) ~ ".bam");
        tmpfiles ~= tmpfilename;
        auto batch_writer = new_bamwriter_from_template(tmpfilename,
                                                        subject_reader, taskpool);
        scope(exit) batch_writer.finish();

        // Write matching reads into tmp
        foreach (read; subject_reader.reads) {
            if (read.name in cache) {
                batch_writer.writeRecord(read);
            }
        }
    }
    writeln("[filter_bam] - combining tmp bam files");
    auto multireader = new MultiBamReader(tmpfiles, taskpool);
    auto writer = new_bamwriter_from_template(outfile, multireader, taskpool);
    scope(exit) writer.finish();

    foreach (read; multireader.reads()) {
        writer.writeRecord(read);
    }

    writeln("[filter_bam] - cleaning up");
    foreach (tmpfile; tmpfiles) {
        remove(tmpfile);
    }
}

int threads = 1;
int MAPQUAL = 30;
int BASEQUAL = 15;
int MINCOV = 10;
auto delete_wdir = true;
string working_dir = ".sidekick.tmp";
string fl_in;  // Input bam file
string hm_out; // Half-mapped out file
string hu_out; // Half-unmapped out file
string au_out; // All-unmapped out file

void main(string[] argv)
{
    // Argument handling
    auto args = getopt(
        argv,
        "threads|t", &threads,
        "mapqual|q", &MAPQUAL,
        "basequal|b", &BASEQUAL,
        "coverage|c", &MINCOV,
        "working_dir|w", &working_dir,
        "delete|d", &delete_wdir,
        std.getopt.config.required, "input|i", &fl_in,
        std.getopt.config.required, "mapped|m", &hm_out,
        std.getopt.config.required, "unmapped|u", &hu_out,
        std.getopt.config.required, "all|a", &au_out);

    writefln("threads = %d", threads);
    writefln("MAPQUAL = %d", MAPQUAL);
    writefln("BASEQUAL = %d", BASEQUAL);
    writefln("MINCOV = %d", MINCOV);

    if (args.helpWanted)
    {
        defaultGetoptPrinter("", args.options);
        import core.stdc.stdlib: exit;
        exit(0);
    }

    // Set up temporary working dir
    mkdirRecurse(working_dir);
    assert(exists(working_dir));

    scope(exit) {
        if (delete_wdir) {
            rmdirRecurse(working_dir);
        }
    }

    // Define paths to tempfiles
    string tmp_mapped = buildPath(working_dir, "tmp_mapped.bam");  // Mapped half of mapped-unmapped pairs go here
    string tmp_mapped_filtered = buildPath(working_dir, "tmp_mapped_filtered.bam");  // Filtered version of above
    string tmp_unmapped = buildPath(working_dir, "tmp_unmapped.bam");  // Unmapped half of mapped-unmapped pairs
    string tmp_unmapped_filtered = buildPath(working_dir, "tmp_unmapped_filtered.bam");
    string tmp_both_1 = buildPath(working_dir, "tmp_both_1.bam");  // Both unmapped, first in pair
    string tmp_both_2 = buildPath(working_dir, "tmp_both_2.bam");  // Both unmapped, second in pair

    // Set up parallelism for bam reader
    auto pool = new TaskPool(threads);
    scope(exit) pool.finish();

    // Open the input bam files with a parallel pool
    auto reader = new BamReader(fl_in, pool);

    // Declare regionlist in this scope
    BamRegion[] regions;

    // 1: Select read pairs where at least one read is unmapped
    // 2: Apply quality filters,
    //    i.e. mapped reads >=MAPQUAL mapping quality, unmapped >=BASEQUAL avg base quality
    // 3: Divert half-mapped, half-unmapped and fully unmapped reads into separate files
    // -caveat: using assoc array 'cache' to keep track of read pairing
    writefln("Searching %s for unmapped reads", fl_in);
    writeln("----");
    {
        auto bar = new shared(ProgressBar)(100000);
        scope(exit) bar.finish();GC.collect();GC.minimize();
        int nreads = 0;

        auto unmapfile = new_bamwriter_from_template(tmp_unmapped, reader, pool);
        scope(exit) unmapfile.finish();

        // auto unmapfile_2 = new_bamwriter_from_template(tmp_unmapped_2, reader, pool);
        // scope(exit) unmapfile_2.finish();

        auto both_first_file = new_bamwriter_from_template(tmp_both_1, reader, pool);
        scope(exit) both_first_file.finish();

        auto both_second_file = new_bamwriter_from_template(tmp_both_2, reader, pool);
        scope(exit) both_second_file.finish();

        auto mapfile = new_bamwriter_from_template(tmp_mapped, reader, pool);
        scope(exit) mapfile.finish();

        auto readGetter = reader.readsWithProgress((lazy float p) { bar.update(p); })
                                .filter!(passes_initial_checks)
                                .filter!(r => passes_quality_checks(r, BASEQUAL, MAPQUAL));

        foreach (read; readGetter) {
            nreads++;
            if (read.is_unmapped) {
                if (read.mate_is_unmapped) {
                    if (read.is_first_of_pair) both_first_file.writeRecord(read);
                    else if (read.is_second_of_pair) both_second_file.writeRecord(read);
                }
                else {
                    unmapfile.writeRecord(read);
                    //writefln("Inserting '%s' into the cache", read.name);
                    //cache.insert(read.name);
                    //assert(cache.contains(read.name));
                }
            }
            else {
                mapfile.writeRecord(read);
            }
        }
    }
    
    // 1: Filter tmp_mapped to only contain reads with a partner in tmp_unmapped => tmp_mapped_filtered
    // 2: Filter tmp_unmapped for reads in tmp_mapped_filtered => tmp_unmapped_filtered
    writeln("Filtering tmp bam files for matching pairs");
    {
        scope(exit) GC.collect();GC.minimize();
        filter_bam(tmp_unmapped, tmp_mapped, pool, working_dir, tmp_mapped_filtered);
        filter_bam(tmp_mapped_filtered, tmp_unmapped, pool, working_dir, tmp_unmapped_filtered);
    }

    // 1: Reopen half-mapped reads file
    // 2: Filter for reads that have a quality passing unmapped pair (in cache)
    //    and divert to filtered file
    // 3: Pileup half-mapped reads and identify regions where depth > MINCOV
    writeln("Checking coverage of filtered reads");
    writeln("----");
    {
        auto bar = new shared(ProgressBar)(10000);
        scope(exit) bar.finish();GC.collect();GC.minimize();

        auto mReader = new BamReader(tmp_mapped_filtered, pool);

        //auto filteredmapfile = new_bamwriter_from_template(tmp_mapped_filtered, mReader, pool);
        //scope(exit) filteredmapfile.finish();

        auto mReadGetter = mReader.readsWithProgress((lazy float p) { bar.update(p); });        

        if (mReadGetter.empty) {
            import std.datetime: Clock;
            writefln("\n[%s] Error - No qualifying reads are available", Clock.currTime.toSimpleString());
            goto END;
        }

        // Now do pileup and identify regions
        foreach (chunk; pileupChunks(mReadGetter, true)) {

            int refid = 0; // Region delimiters
            ulong left = 0;
            ulong right = 0;
            bool initialised = false;

            foreach (column; chunk) {
                if (column.coverage >= MINCOV) {
                    // writefln("%d",
                        // column.coverage);
                    // foreach (read; column.reads) {
                        // writefln("%s\t%d\t%d",
                            // read.name, read.position, read.position + read.basesCovered());
                    // }
                    if (!initialised) {
                        //writeln("init...");
                        initialised = true;
                        refid = column.ref_id;
                        left = column.position;
                        right = column.position;
                    }
                    right = column.position;
                }
                else if (initialised) {
                    auto region = BamRegion(refid, to!uint(left), to!uint(right));
                    regions ~= region;
                    initialised = false;
                    refid = left = right = 0;
                }
            }
            // Reached end of chunk without finalising, so do it here
            if (initialised) {
                auto region = BamRegion(refid, to!uint(left), to!uint(right));
                regions ~= region;
                initialised = false;
                refid = left = right = 0;
            }
        }
    }
    // writeln(regions);

    if (regions.length == 0) {
        import std.datetime: Clock;
        writefln("\n[%s] Error - No qualifying regions were found.", Clock.currTime.toSimpleString());
        goto END;
    }

    // 1: Reopen cached unmapped reads
    // 2: Filter for reads paired to coverage-passing mapped reads
    writeln("Writing half-unmapped reads tied to high coverage areas");
    {
        reader = new BamReader(tmp_unmapped_filtered, pool);
        int i=0;
        auto region = regions[i];

        auto writer = new_bamwriter_from_template(hu_out, reader, pool);
        scope(exit) writer.finish();GC.collect();GC.minimize();

        foreach (read; reader.reads) {
            if (overlaps(read, region, true)) {
                writer.writeRecord(read);
            }
            else if (region.fullyLeftOf(read.ref_id, read.mate_position)) {
                i++;
                if (i >= regions.length) break;
                region = regions[i];
            }
        }
    }

    // 1: Reopen cached mapped reads
    // 2: Filter for coverage-passing reads
    writeln("Writing half-mapped reads tied to high coverage areas");
    {
        reader = new BamReader(tmp_mapped_filtered, pool);
        int i=0;

        auto region = regions[i];

        auto writer = new_bamwriter_from_template(hm_out, reader, pool);
        scope(exit) writer.finish();GC.collect();GC.minimize();

        foreach (read; reader.reads) {
            if (overlaps(read, region, false)) {
                writer.writeRecord(read);
            }
            else if (region.fullyLeftOf(read.ref_id, read.mate_position)) {
                i++;
                if (i >= regions.length) break;
                region = regions[i];
            }
        }
    }

    writefln("Completed. Half mapped reads written to %s.", hm_out);
    writefln("           Half unmapped reads written to %s.", hu_out);
    writefln("           Fully unmapped reads written to %s.", au_out);
END:
    import core.thread;
    Thread.sleep( dur!("seconds")( 30 ) );
    writeln("Done.");
}
