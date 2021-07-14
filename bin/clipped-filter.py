#!/usr/bin/env python

from collections import defaultdict
import itertools
import pysam

# Requires input BAM stream, primary alns only
# Note this doesn't have to be name-sorted but it's a lot easier
# on memory if it is :)

# logic:

# if either R1 or R2 is clipped (has a 3' or 5' soft clip) keep both reads

# Basic code template came from: https://www.biostars.org/p/306041/
def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: None)
	
	# capture all primary alns with same read ID (R1 and R2) in dicts
    for read in bam:
        qname = read.query_name
        if qname not in read_dict:
            read_dict[qname] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname]
            else:
                yield read_dict[qname], read
            del read_dict[qname]

def is_clipped(read, cutoff):
	"""
	Test if a read is soft-clipped but also lacks supplementary alignments (not chimeric)
	"""
	return not(read.has_tag('SA')) and read.get_cigar_stats()[0][4] >= cutoff

save = pysam.set_verbosity(0)
bam = pysam.AlignmentFile('-', 'rb', threads = 4, index_filename=None)
pysam.set_verbosity(save)

outfile = pysam.AlignmentFile("-", "wb", template=bam)

for read1, read2 in read_pair_generator(bam):
	if is_clipped(read1, 20) or is_clipped(read2, 20):
		outfile.write(read1)
		outfile.write(read2)
