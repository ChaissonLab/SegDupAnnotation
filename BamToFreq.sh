#!/usr/bin/env bash
rgn=`echo $2 | tr ":-" "__"`


bamToFreq $1 $2 $3 --covbed | \
 awk '{ if ($3 > 5) print;}' | \
 awk '{ print $1"\t"$2"\t"$2+1"\t"$3;}' | \
 bedtools merge -c 4,4,4 -o mean,stdev,mode > $4.$rgn.rgn 2> /dev/null

true
