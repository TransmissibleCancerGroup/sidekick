import bio.bam.read: compareCoordinates;
import bio.bam.region: BamRegion;
import bio.bam.reader: BamReader;
import bio.bam.writer: BamWriter;
import bio.bam.pileup: pileupChunks;
import utils.progressbar: ProgressBar;
import std.algorithm: filter, reduce, sort;
import std.conv: to;
import std.file: remove;
import std.functional: toDelegate;
import std.parallelism: TaskPool;
import std.range: iota, tee;
import std.stdio: writeln, writefln;

bool passes_initial_checks(R)(const R r) {
	return (!r.is_duplicate &&
		    r.is_paired &&
			!r.proper_pair &&
			!r.failed_quality_control &&
			!r.is_supplementary &&
			!r.is_secondary_alignment &&
			(r.is_unmapped || r.mate_is_unmapped));
}

// Use with std.functional.partial to fix the qual values
bool passes_quality_checks(R)(const R r, int base_qual, int map_qual) {
	return r.is_unmapped ? avg_base_quality(r) >= base_qual : r.mapping_quality >= map_qual;
}

double avg_base_quality(R)(const ref R r) {
	auto qual = reduce!"a + b"(0.0, r.base_qualities);
	return qual / r.base_qualities.length;
}

void write_unmapped(Read, FileHandle)(ref Read r, ref FileHandle file) {
	if (r.is_unmapped) file.writeRecord(r);
}

bool overlaps(Read, Region)(Read read, Region region, bool use_mate_info=false) {
	auto read_start = use_mate_info ? read.mate_position : read.position;
	auto read_ref_id = use_mate_info ? read.mate_ref_id : read.ref_id;
	return (read_ref_id == region.ref_id && 
	    	read_start <= region.end && 
			read_start + (use_mate_info ? read.sequence_length : read.basesCovered()) > region.start);
}

auto new_bamwriter_from_template(string filename, ref BamReader reader, ref TaskPool taskpool) {
	// Remember to call scope(exit) writer.finish() on the returned writer
	auto writer = new BamWriter(filename, -1, taskpool);
	writer.writeSamHeader(reader.header);
	writer.writeReferenceSequenceInfo(reader.reference_sequences);
	return writer;
}

void main(string[] argv)
{
	// Argument handling (TODO: do properly later w/ getopt)
	if (argv.length != 6)
	{
		writefln("Usage: %s <threads> <bamfile:in> <half_mapped_reads_bamfile:out> <half_unmapped_reads_bamfile:out> <all_unmapped_reads_bamfile:out>", argv[0]);
		import std.c.process: exit;
		exit(-1);
	} 
	int nthreads = to!int(argv[1]);
	string fl_in = argv[2];
	string hm_out = argv[3];
	string hu_out = argv[4];
	string au_out = argv[5];

	// Collect thresholds (TODO: make getopt options)
	int MAPQUAL = 30;
	int BASEQUAL = 10;
	int MINCOV = 1;

	// Set up parallelism for bam reader
	auto pool = new TaskPool(nthreads);
	scope(exit) pool.finish();

	// Open the input bam files with a parallel pool
	auto reader = new BamReader(fl_in, pool);
	
	// Declare regionlist in this scope
	BamRegion[] regions;
	bool[string] cache;	

	// 1: Select read pairs where at least one read is unmapped
	// 2: Apply quality filters, 
	//    i.e. mapped reads >=MAPQUAL mapping quality, unmapped >=BASEQUAL avg base quality
	// 3: Divert half-mapped, half-unmapped and fully unmapped reads into separate files 
	// -caveat: using assoc array 'cache' to keep track of read pairing
	{
		writefln("Searching %s for unmapped reads", fl_in);
		auto bar = new shared(ProgressBar)(100000);
		scope(exit) bar.finish();
		int nreads = 0;
		auto unmapfile = new_bamwriter_from_template("tmp_unmapped.bam", reader, pool);
		scope(exit) unmapfile.finish();
		auto doubleunmapfile = new_bamwriter_from_template(au_out, reader, pool);
		scope(exit) doubleunmapfile.finish();
		auto mapfile = new_bamwriter_from_template("tmp_mapped.bam", reader, pool);
		scope(exit) mapfile.finish();
		
		auto readGetter = reader.readsWithProgress((lazy float p) { bar.update(p); })
								.filter!(passes_initial_checks)
								.filter!(r => passes_quality_checks(r, BASEQUAL, MAPQUAL));
	
		foreach (read; readGetter) {
			nreads++;
			if (read.is_unmapped) {
				if (read.mate_is_unmapped) {
					doubleunmapfile.writeRecord(read);
				}
				else {
					unmapfile.writeRecord(read);
					cache[read.name] = true;
				}
			}
			else {
			    mapfile.writeRecord(read);
			}
		}
	}

	// 1: Reopen half-mapped reads file
	// 2: Filter for reads that have a quality passing unmapped pair (in cache)
    //    and divert to filtered file
	// 3: Pileup half-mapped reads and identify regions where depth > MINCOV
	{
		writeln("Checking coverage of filtered reads");
		auto bar = new shared(ProgressBar)(1000);
		scope(exit) bar.finish();
		
		auto mReader = new BamReader("tmp_mapped.bam", pool);
		scope(exit) {
			remove("tmp_mapped.bam");
		    remove("tmp_mapped.bam.bai");
		}
		
		auto filteredmapfile = new_bamwriter_from_template("tmp_mapped_filtered.bam", mReader, pool);
		scope(exit) filteredmapfile.finish();
		
		auto mReadGetter = mReader.readsWithProgress((lazy float p) { bar.update(p); })
							.filter!(r => r.name in cache)
							.tee!(r => filteredmapfile.writeRecord(r));
		
		// Now do pileup and identify regions
		foreach(chunk; pileupChunks(mReadGetter, true)) {
			
			int refid = 0; // Region delimiters
			ulong left = 0;
			ulong right = 0;
			bool initialised = false;
			
			foreach(column; chunk) {
				if (column.coverage > MINCOV) {
					if (!initialised) {
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
		}
		cache = cache.init; // kill the cache
	}

	// 1: Reopen cached unmapped reads
	// 2: Filter for reads paired to coverage-passing mapped reads
	{
		writeln("Writing half-unmapped reads tied to high coverage areas");
		reader = new BamReader("tmp_unmapped.bam", pool);
		int i=0;
		auto region = regions[i];
		
		auto writer = new_bamwriter_from_template(hu_out, reader, pool);
		scope(exit) {
			remove("tmp_unmapped.bam");
			remove("tmp_unmapped.bam.bai");
			writer.finish();
		}

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
	{
		writeln("Writing half-mapped reads tied to high coverage areas");
		reader = new BamReader("tmp_mapped_filtered.bam", pool);
		int i=0;
		auto region = regions[i];
		
		auto writer = new_bamwriter_from_template(hm_out, reader, pool);
		scope(exit) {
			remove("tmp_mapped_filtered.bam");
			remove("tmp_mapped_filtered.bam.bai");
			writer.finish();
		}
		
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
}
