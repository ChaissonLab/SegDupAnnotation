#!/usr/bin/env bash
resDups=$1
bamFile=$2

while read line; do
		rgn1=`echo $line | awk '{ print $1":"$2"-"$3;}'`
    d1=`samtools depth $bamFile -r $rgn1 | cut -f 3 | stats | cut -f 1`
		rgn2=`echo $line | awk '{ print $4":"$5"-"$6;}'`
    d2=`samtools depth $bamFile -r $rgn2 | cut -f 3 | stats | cut -f 1`
		l1=`echo $line | awk '{ print $3-$2;}'`		
		l2=`echo $line | awk '{ print $6-$5;}'`		
		echo -e $rgn1"\t"$rgn2"\t"$d1"\t"$d2"\t"$l1"\t"$l2
done < $resDups
