#!/usr/bin/env bash

table=$1
assembly=$2
genes=$3
bam=$4;
output=$5;

samtools view -H $bam > $output

while read line; do 

gene=`echo $line | tr " " "\t" | cut -f 4`;
orig=`echo $line | awk '{ print $13":"$14"-"$15;}'`
dup=`echo $line | awk '{ print $16":"$17"-"$18;}'`


mkdir -p realign

samtools faidx $genes "$gene" > realign/gene.fasta

samtools faidx $assembly $orig > realign/orig.fasta

minimap2 -a -x splice realign/orig.fasta realign/gene.fasta  > realign/aln.sam

ch=`echo $line | tr " " "\t" | cut -f 13`
pos=`echo $line | tr " " "\t" |  cut -f 14`
cat realign/aln.sam | awk  -vchr="$ch" -vpos="$pos" '{ if (substr($1,0,1) != "@") { $3 = chr; $4+=pos-1; print; } }' | tr " " "\t" >> $output


samtools faidx $assembly $dup > realign/ref.fasta

minimap2 -a -x splice realign/ref.fasta realign/gene.fasta  > realign/aln.sam

ch=`echo $line | tr " " "\t" | cut -f 16`
pos=`echo $line | tr " " "\t" |  cut -f 17`
echo $ch
echo $pos

cat realign/aln.sam | awk  -vchr="$ch" -vpos="$pos" '{ if (substr($1,0,1) != "@") { $3 = chr; $4+=pos-1; print; } }' | tr " " "\t" >> $output

done < $table
