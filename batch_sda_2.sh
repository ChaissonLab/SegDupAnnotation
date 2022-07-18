#!/bin/bash

# Keon Rabbani
# krabbani@usc.edu
# 12/28/2021
# Updated: 05/24/2022

# Purpose: Run SegDupAnalysis on genomes given VGP's species name and code.

set -euo pipefail

numArgsRequired=1

if [ "$#" -lt "$numArgsRequired" ] # is there an argument?
then
    echo "error: too few arguments, you provided $#, $numArgsRequired required"
    echo "usage: $0 speciesCode.txt"
    echo "speciesCode.txt has at least 2 columns: col1=speciesName, col2=speciesCode"
    exit 1
fi

echo "------ running $0 ------"

# Save Time
date
res1=$(date +%s.%N)

# Load CARC Modules
#module load usc gcc aws-cli boost

# set & save vars
export speciesCodenameFile="$1"

export dataDirPath="/scratch1/krabbani/data/"
export analysisDirPath="/scratch2/krabbani/analysis_finch/"
export tmpDirPath="/scratch2/krabbani/sda_tmp/"
export jsonTemplatePath="/scratch2/krabbani/analysis/sd_analysis_template.json"
export initialDir=$(pwd)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib

export email="krabbani@usc.edu"

# Loop through each line of input file
while IFS=$'\t' read -r speciesName speciesCode remainder
do
    echo "---- starting: $speciesName ($speciesCode) ----"
    # Create Dirs
    echo "-- creating dirs --"
    mkdir -p "$dataDirPath"/raw_"$speciesCode"
    mkdir -p "$analysisDirPath"/sda_"$speciesCode"
    mkdir -p "$tmpDirPath"/tmp_"$speciesCode"_tmp # TODO FINCH

#ssh krabbani@hpc-transfer1.usc.edu /bin/bash <<EndOfRemoteExec
ssh hpc-transfer1 /bin/bash <<EndOfRemoteExec
    ## Data Download moved to another script to be run on hpc-transfer.usc.edu
    # Download Data
    echo "-- downloading data --"
    cd "$dataDirPath"/raw_"$speciesCode"
    aws s3 --no-sign-request sync s3://genomeark/species/"$speciesName"/"$speciesCode"/genomic_data/pacbio/ . --exclude "*" --include "*.bam"
    aws s3 --no-sign-request sync s3://genomeark/species/"$speciesName"/"$speciesCode"/assembly_curated/ . --exclude "*" --include "*.pri.cur.*.fasta*"
EndOfRemoteExec

    # Determine correct assembly file to use
    echo "-- picking assembly to use --"
    echo "- NOTE Verify correct assembly file chosen!"
    echo "- Method chosen works in most but not all cases."
    cd "$dataDirPath"/raw_"$speciesCode"
    #export assemblyFileGz=$(ls | grep fasta | grep ".pri.cur.*.fasta" | sort | tail -1)
    export assemblyFileGz=$( ls | grep fasta | grep ".pri.cur.*.fasta" | grep -v "gfastats.tsv" | sort -t '.' -k4 -n | tail -1)
    if [ ${assemblyFileGz:(-3):3} == ".gz" ] && [ ! -f ${assemblyFileGz/%.gz/} ] # if assemblyFileGz ends with gz && an unzipped version does not exist
    then
        echo "Unzipping assembly"
        gunzip "$assemblyFileGz"
        export assemblyFile=${assemblyFileGz%.gz}
    else
        export assemblyFile="$assemblyFileGz"
    fi
    echo "Assembly file chosen: $assemblyFile"

    # Create json config file
    echo "-- creating json --"
    cd "$analysisDirPath"/sda_"$speciesCode"
    cp "$jsonTemplatePath" "$analysisDirPath"/sda_"$speciesCode"/sd_analysis.json
    sed -i "s|##SPECIES_CODE##|$speciesCode|" sd_analysis.json
    sed -i "s|##ASSEMBLY##|$assemblyFile|" sd_analysis.json
    ls "$dataDirPath"/raw_"$speciesCode"/ | grep ".bam$" | sed -e "s|^|       \"/scratch1/krabbani/data/raw_$speciesCode/|" | sed -e 's|$|",|' > formattedBamFileNames.txt
    sed -i 3rformattedBamFileNames.txt "$analysisDirPath"/sda_"$speciesCode"/sd_analysis.json # add contents of 'formattedBamFileNames.txt' after 3rd line of json file
    rm formattedBamFileNames.txt

    # Run snakemake
    cd "$analysisDirPath"/sda_"$speciesCode"
    snakemake -p -s /project/mchaisso_100/cmb-16/krabbani/git/SegDupAnnotation/SegDupAnalysis.snakefile -c 16 --cluster "{params.grid_opts}" -j 300

done < "$speciesCodenameFile"


function finish {
    # Clean Files
    #rm -rf "$tmpDirPath"/tmp_"$speciesCode"/*

    echo

    # Print Time
    res2=$(date +%s.%N) # record script end time
    dt=$(echo "$res2 - $res1" | bc)
    dd=$(echo "$dt/86400" | bc)
    dt2=$(echo "$dt-86400*$dd" | bc)
    dh=$(echo "$dt2/3600" | bc)
    dt3=$(echo "$dt2-3600*$dh" | bc)
    dm=$(echo "$dt3/60" | bc)
    ds=$(echo "$dt3-60*$dm" | bc)
    date
    printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds

    # Send Email
    printf "$0 complete.\nTotal runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds | mail -s "[CARC-POST] $0 done" "$email"

    cd "$initialDir"
}
trap finish EXIT INT QUIT TERM
