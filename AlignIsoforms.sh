#!/usr/bin/env bash
bedFile=$1
gencode=$2
genome=$3
#set -v
#set -x
mkdir -p testing
for gene in `cut -f 4 $bedFile | uniq`; do 
		res=$gene
		for rgn in `grep "$gene" $bedFile | awk '{ print $13"\t"$14"\t"$15; print $16"\t"$17"\t"$18;}' | bedtools sort | bedtools merge  | awk '{ print $1":"$2"-"$3;}'`; do

				samtools faidx $gencode $gene > testing/transcript.fa
				samtools faidx $genome $rgn > testing/genome.fa
				dv=`minimap2 -x splice -O 8,20 testing/genome.fa testing/transcript.fa | awk '{ for (i=1; i <= NF; ++i) { if (substr($i,0,2) == "dv") print(substr($i,6)) } }'`
				res="$res\t$rgn\t $dv"
				done
		echo -e $res
		done

		
