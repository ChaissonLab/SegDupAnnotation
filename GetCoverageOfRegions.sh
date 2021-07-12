#!/usr/bin/env bash
bedFile=$1
covFile=$2
meanCov=`cat $3`


while read line; do
    rgn=`echo $line | awk '{ print $1":"$2"-"$3;}'`
    tabix $covFile  $rgn | cut -f 4 | stats | awk -v m=$meanCov '{ print $1/m;}'
done  < $bedFile > $bedFile.txt
   
