#!/usr/bin/env bash
bed=$1
assembly=$2
mkdir -p dynamic_mask
cut -f 1-6 $bed | awk '{ if ($3-$2 > 200 && $3-$2 < 400 && $6-$5 > 200 && $6-$5 <400) print; }' | bedtools merge -c 1 -o count | awk '{ if ($4 > 20) print;}'b > dynamic_mask/small.bed & 

cut -f 1-6 $bed | awk '{ if ($3-$2 > 400 && $3-$2 < 800 && $6-$5 > 400 && $6-$5 < 800 ) print; }' | bedtools merge -c 1 -o count | awk '{ if ($4 > 20) print;}' > dynamic_mask/medium.bed &

cut -f 1-6 $bed | awk '{ if ($3-$2 > 800 && $3-$2 < 1200 && $6-$5 > 800 && $6-$5 < 1200) print; }' | bedtools merge -c 1 -o count | awk '{ if ($4 > 20) print;}' > dynamic_mask/large.bed &

cut -f 1-6 $bed | awk '{ if ($3-$2 > 1200 && $3-$2 < 5000 && $6-$5 > 1200 && $6-$5 < 5000) print; }' | bedtools merge -c 1 -o count | awk '{ if ($4 > 20) print;}' > dynamic_mask/very_large.bed &

cut -f 1- 6 $bed | awk '{ if ($3-$2 > 5000 && $3-$2 < 8000 && $6-$5 > 5000 && $6-$5 < 8000) print; }' | bedtools merge -c 1 -o count | awk '{ if ($4 > 20) print;}' > dynamic_mask/LINE.bed &

wait

#bedtools intersect -a $bed -b dynamic_mask/small.bed -wa -f 0.7 -sorted | awk '{ print $1"\t"$2"\t"$3; print $4"\t"$5"\t"$6;}' | bedtools sort | bedtools merge > dynamic_mask/all_small.bed &

bedtools intersect -a $bed -b dynamic_mask/medium.bed -wa -f 0.7 -sorted | awk '{ print $1"\t"$2"\t"$3; print $4"\t"$5"\t"$6;}' | bedtools sort | bedtools merge > dynamic_mask/all_medium.bed &

bedtools intersect -a $bed -b dynamic_mask/large.bed -wa -f 0.7 -sorted | awk '{ print $1"\t"$2"\t"$3; print $4"\t"$5"\t"$6;}' | bedtools sort | bedtools merge > dynamic_mask/all_large.bed &

bedtools intersect -a $bed -b dynamic_mask/very_large.bed -wa -f 0.7 -sorted | awk '{ print $1"\t"$2"\t"$3; print $4"\t"$5"\t"$6;}' | bedtools sort | bedtools merge > dynamic_mask/all_very_large.bed &

bedtools intersect -a $bed -b dynamic_mask/LINE.bed -wa -f 0.7 -sorted | awk '{ print $1"\t"$2"\t"$3; print $4"\t"$5"\t"$6;}' | bedtools sort | bedtools merge > dynamic_mask/all_LINE.bed &
wait

