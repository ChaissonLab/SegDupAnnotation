#!/bin/bash



sum=/home/cmb-16/mjc/rdagnew/summerproj
bam=$1
fai=$2



#intersectBed -c -a <(bedtools makewindows -b <( awk 'BEGIN{OFS="\t"}{print $1,0,$2}' $fai) -w 100) -b samBed.$bam.bed.f > coverage.$bam.bed

awk 'BEGIN{OFS="\t";c=0;sum=0;} sum=sum+$4;c=c+1;END{print sum/c;}' coverage.$bam.bed |tail -1> mean.$bam.txt

mean=$(cat mean.$bam.txt)
echo $mean

#viterbi2.cpp
for r in ` cat $fai|cut -f 1`;do echo $r; $sum/hmm/viterbi <(grep -w $r coverage.$bam.bed |cut -f 4 ) $mean $bam.$r; done

for r in ` cat $fai|cut -f 1`;do echo $r; cat $bam.$r.viterout.txt >> $bam.viterout.txt;done

paste <(cat coverage.$bam.bed) <(cat $bam.viterout.txt) > $bam.viterout.bed

#plot every 50000 points ~ 5MB per plot
/home/cmb-16/mjc/rdagnew/anaconda3/bin/Rscript $sum/plot.HMM.noclip.R $bam.viterout.bed $bam 50000 mean.$bam.txt
