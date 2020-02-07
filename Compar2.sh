#!/usr/bin/env bash

read -r line
echo $line
ra=`echo $line | awk '{ print $1":"$2"-"$3;}'`
rb=`echo $line | awk '{ print $4":"$5"-"$6;}'`
echo $line > testing/rep.bed
samtools faidx assembly.repeat_masked.sd.fasta $ra > testing/a.fasta
samtools faidx assembly.repeat_masked.sd.fasta $rb > testing/b.fasta

minimap2 --cs=long testing/a.fasta testing/b.fasta | paftools.js view  -

