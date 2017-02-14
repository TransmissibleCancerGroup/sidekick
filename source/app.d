import bio.bam.reader;
import bio.bam.pileup;
import std.algorithm: filter, reduce;
import std.conv;
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

// void progress(lazy float p) {
// 	static uint n;
// 	if(++n % 100 == 0) {
// 		writeln(p)
// 	}
// }

void main(string[] argv)
{
	// Argument handling (do properly later)
	if (argv.length != 3)
	{
		writefln("Usage: %s <threads> <bamfile>", argv[0]);
		import std.c.process: exit;
		exit(-1);
	} 
	int nthreads = to!int(argv[1]);
	string filename = argv[2];

	// Set up parallelism for bam reader
	auto pool = new TaskPool(nthreads);
	scope(exit) pool.finish();

	// Open the bam file with a parallel reader
	auto bam = new BamReader(filename, pool);

	// Iterate over the bam file
	BamRead[] both_unmapped;
	BamRead[] mate_mapped;
	BamRead[] i_am_mapped_my_mate_is_not;
	BamRead[] [string] data;

	foreach (read; filter!(passes_initial_checks)(bam.reads)) {
		data[read.name] ~= read;
		if (read.is_unmapped && avg_base_quality(read) > 15) {
			if (read.mate_is_unmapped) {
			    both_unmapped ~= read;
		    }
			else {
				mate_mapped ~= read;
			}
		}
		else if (read.mate_is_unmapped && read.mapping_quality > 30) {
			i_am_mapped_my_mate_is_not ~= read;
		}
	}
	
	writefln("Found %d half-unmapped pairs, and %d fully unmapped pairs",
	    mate_mapped.length, both_unmapped.length);
	writefln("Found %d mapped partners of half-unmapped pairs",
	    i_am_mapped_my_mate_is_not.length);
	writefln("Found %d total", data.length);

	foreach (key; data.keys) {
		writefln("%s\t%d", key, data[key].length);
	}
	// TODO: Collect the mapped mates of our mate_mapped reads and build pileup
	auto pileup = makePileup(i_am_mapped_my_mate_is_not);
	
	foreach (column; pileup) {
		if (column.reference_base == 'N') continue;
		writeln("Column position: ", column.position);
		writeln("    Ref.base: ", column.reference_base); // extracted from MD tags
		writeln("    Coverage: ", column.coverage);
	}

	writeln("Done.");
}
