import bio.bam.read: compareCoordinates;
import bio.bam.reader;
import bio.bam.writer;
import bio.bam.pileup;
import utils.progressbar;
import std.algorithm: filter, reduce, sort;
import std.conv;
import std.functional: toDelegate;
import std.parallelism;
import std.stdio;

bool passes_initial_checks(R)(const ref R r) {
	return (!r.is_duplicate &&
		    r.is_paired &&
			!r.proper_pair &&
			!r.failed_quality_control &&
			!r.is_supplementary &&
			!r.is_secondary_alignment &&
			(r.is_unmapped || r.mate_is_unmapped));
}

double avg_base_quality(R)(const ref R r) {
	auto qual = reduce!"a + b"(0.0, r.base_qualities);
	return qual / r.base_qualities.length;
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
	int MAPQUAL = 50;
	int BASEQUAL = 15;
	int MINCOV = 1;

	// Set up parallelism for bam reader
	auto pool = new TaskPool(nthreads);
	scope(exit) pool.finish();

	// Open the input bam files with a parallel pool
	auto reader = new BamReader(fl_in, pool);

	// Iterate over the bam file
	BamRead[] [string] data;

	auto bar = new shared(ProgressBar)(100000);
	foreach (read; (reader.readsWithProgress((lazy float p) { bar.update(p); }))) {
		if (passes_initial_checks(read)) data[read.name] ~= read;
	}
	bar.finish();
	
	writefln("Found %d total", data.length);

	foreach (key; data.keys) {
		sort!("a.is_first_of_pair")(data[key]);
		if (data[key].length > 1) assert(data[key][0].is_first_of_pair && data[key][1].is_second_of_pair);
	}

	// Get reads for pileup
	BamRead[] mapped;

	foreach(readList; data.values) {
		if (readList.length < 2) continue;
		if (readList[1].is_unmapped && readList[0].mapping_quality > MAPQUAL && avg_base_quality(readList[1]) > BASEQUAL) {
			mapped ~= readList[0];
		}
		else if (readList[0].is_unmapped && readList[1].mapping_quality > MAPQUAL && avg_base_quality(readList[0]) > BASEQUAL) {
			mapped ~= readList[1];
		}
	}
	sort!compareCoordinates(mapped);
	writefln("Found %d reads for pileup.", mapped.length);
	
	// Pileup to select reads that are in regions that achieve a certain coverage
	if (mapped.length > 0) {
		
		// Now open the output bam, reusing the Task Pool
		auto half_mapped_writer = new BamWriter(hm_out, -1, pool);
		scope(exit) half_mapped_writer.finish();
		half_mapped_writer.writeSamHeader(reader.header);
		half_mapped_writer.writeReferenceSequenceInfo(reader.reference_sequences);

		auto half_unmapped_writer = new BamWriter(hu_out, -1, pool);
		scope(exit) half_unmapped_writer.finish();
		half_unmapped_writer.writeSamHeader(reader.header);
		half_unmapped_writer.writeReferenceSequenceInfo(reader.reference_sequences);

		auto all_unmapped_writer = new BamWriter(au_out, -1, pool);
		scope(exit) all_unmapped_writer.finish();
		all_unmapped_writer.writeSamHeader(reader.header);
		all_unmapped_writer.writeReferenceSequenceInfo(reader.reference_sequences);
		
		foreach(chunk; pileupChunks(mapped, true)) {
			foreach (column; chunk) {
				if (column.coverage > MINCOV) {
					foreach (read; column.reads) {
						if (read.name in data) {
							BamRead[] readlist = data[read.name];
							if (readlist[0].is_unmapped && readlist[1].is_unmapped) {
								all_unmapped_writer.writeRecord(readlist[0]);
								all_unmapped_writer.writeRecord(readlist[1]);
							}
							else if (!readlist[0].is_unmapped && readlist[1].is_unmapped) {
								half_mapped_writer.writeRecord(readlist[0]);
								half_unmapped_writer.writeRecord(readlist[1]);
							} 
							else if (readlist[0].is_unmapped && !readlist[1].is_unmapped) {
								half_unmapped_writer.writeRecord(readlist[0]);
								half_mapped_writer.writeRecord(readlist[1]);
							}
							data.remove(read.name);
						}
					}
				}
			}
		}
		writefln("Completed. Half mapped reads written to %s.", hm_out);
		writefln("           Half unmapped reads written to %s.", hu_out);
		writefln("           Fully unmapped reads written to %s.", au_out);
	}
	else {
		writeln("Completed, but no unmapped reads found.");
	}
}
