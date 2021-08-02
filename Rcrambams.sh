#!/usr/bin/env bash
samtools view ref_aligned.bam -C -@ 4 -T assembly.hg38.fa -o ref_aligned.cram