#!/usr/bin/env bash
if [[ $1 == *".fastq.gz" ]]; then
		zcat $1
elif [[ $1 == *".bam" ]]; then
		samtools view -hF 2304 $1 | samtools view -h | samtools fastq -
elif [[ $1 == *".fastq" ]]; then
		cat $1
fi

