from __future__ import print_function
# Code to process a bam file and extract reads with left behind
# sidekicks (mapped reads with unmapped mates)

import pysam

def hero(read):
    return (read.is_paired and
            not read.is_unmapped and
            read.mate_is_unmapped and
            not read.is_supplementary and
            read.mapping_quality > 50)

def sidekick(read):
    return (read.is_paired and
            read.is_unmapped and
            not read.mate_is_unmapped and
            not read.is_supplementary)
