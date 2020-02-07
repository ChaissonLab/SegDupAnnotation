#!/usr/bin/env bash

while read line; do
rgna=`echo $line | awk '{ print $1":"$2"-"$3;}'`
rgnb=`echo $line | awk '{ print $4":"$5"-"$6;}'`

samtools faidx assembly.orig.fasta $rgna > testing/a.fasta
samtools faidx assembly.orig.fasta $rgnb > testing/b.fasta

blasrmc testing/a.fasta testing/b.fasta -m 0 
done < $1
