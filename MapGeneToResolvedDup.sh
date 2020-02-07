#!/usr/bin/env bash

while read line; do 
		transcript=`echo $line | cut -f 4`
		dupRegion=`echo $line | awk '{ print 

done < $1



