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

void progress(lazy float p) {
	static uint n;
	if(++n % 100 == 0) {
		writeln(p);
	}
}

void main(string[] argv)
{
	// Argument handling (do properly later)
	if (argv.length != 4)
	{
		writefln("Usage: %s <threads> <bamfile:in> <bamfile:out>", argv[0]);
		import std.c.process: exit;
		exit(-1);
	} 
	int nthreads = to!int(argv[1]);
	string filename = argv[2];

	// Set up parallelism for bam reader
	auto pool = new TaskPool(nthreads);
	scope(exit) pool.finish();

	// Open the input and output bam files with a parallel pool
	auto reader = new BamReader(filename, pool);
	auto writer = new BamWriter(argv[3], -1, pool);
	scope(exit) writer.finish();
	writer.writeSamHeader(reader.header);
	writer.writeReferenceSequenceInfo(reader.reference_sequences);

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
		if (readList[1].is_unmapped && readList[0].mapping_quality > 50 && avg_base_quality(readList[1]) > 15) {
			mapped ~= readList[0];
		}
		else if (readList[0].is_unmapped && readList[1].mapping_quality > 50 && avg_base_quality(readList[0]) > 15) {
			mapped ~= readList[1];
		}
	}
	sort!compareCoordinates(mapped);
	writefln("Found %d reads for pileup.", mapped.length);
	
	// Pileup to select reads that are in regions that achieve a certain coverage
	foreach(chunk; pileupChunks(mapped, true)) {
		foreach (column; chunk) {
			if (column.coverage > 0) {
				foreach (read; column.reads) {
					if (read.name in data) {
						foreach (cached_read; data[read.name]) {
							writer.writeRecord(cached_read);
						}
						data.remove(read.name);
					}
				}
			}
		}
	}
	writeln("Done.");
}
