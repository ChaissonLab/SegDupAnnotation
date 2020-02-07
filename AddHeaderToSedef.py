#!/usr/bin/env python
import sys
print("#"+"\t".join(["chr1", "start1", "end1","chr2", "start2","end2","name","score","strand1","strand2","max_len","aln_len", "comment","aln_len","indel_a","indel_b","alnB","matchB","mismatchB","transitionsB","transversions","fracMatch","fracMatchIndel","jck","k2K","aln_gaps","uppercaseA","uppercaseB","uppercaseMatches","aln_matches","aln_mismatches","aln_gaps","aln_gap_bases","cigar","filter_score"]))
for line in sys.stdin:
    sys.stdout.write(line)
