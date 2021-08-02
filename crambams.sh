#!/usr/bin/env bash
samtools view assembly.bam -C -@ 4 -T assembly.orig.fasta -o assembly.cram