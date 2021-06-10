import os
import tempfile
import subprocess
import os.path
import json

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)


#config("hmm_caller.json")
configfile: "sd_analysis.json"

outdir="hmm"

fai= open(config["asm"]+".fai")
contigs = [l.split()[0].strip().replace("|","_") for l in fai]
bamt = config["bam"]
bam = bamt.split("/")[-1]
prefix_bam = os.path.splitext(bam)[0]

rule all:
    input:
        bed="hmm/cov.bed",
        no_subreadbed="hmm/cov.no_subread.bed",
        bins="hmm/coverage.bins.bed.gz",
        avg="hmm/mean_cov.txt",
        vitterout=expand("hmm/{ctg}.viterout.txt", ctg=contigs),
        cn=expand("hmm/copy_number.{ctg}.bed", ctg=contigs),
        allCN="hmm/copy_number.tsv",
#        plot=expand("hmm/{bm}.noclip.pdf",bm=prefix_bam),

rule MakeCovBed:
    input:
        bam=config["bam"],
    output:
        bed=protected("hmm/cov.bed"),
    params:
        sd=SD,
    shell:"""
mkdir -p hmm
samtools view -q 10 -F 2304 -@ 3 {input.bam} | {params.sd}/hmm/samToBed /dev/stdin/ --useH --flag   > {output.bed}
"""

if config['index_params']==" -CLR":    
    rule FilterSubreads:
        input:
            bed="hmm/cov.bed",
        output:
            covbed="hmm/cov.no_subread.bed",
        params:
            sd=SD,
        shell:"""
    {params.sd}/RemoveRedundantSubreads.py --input {input.bed} |sort -k1,1 -k2,2n > {output.covbed}
    """
else:
    rule FilterSubreads0:
        input:
            bed="hmm/cov.bed",
        output:
            covbed="hmm/cov.no_subread.bed",
        params:
            sd=SD,
        shell:"""
cd hmm; ln -s cov.bed cov.no_subread.bed
    """
rule MakeIntersect:
    input:
        bed="hmm/cov.no_subread.bed",
    output:
        bins=protected("hmm/coverage.bins.bed.gz"),
    params:
        asm=config["asm"],
        sd=SD
    shell:"""
{params.sd}/BedToCoverage.py {input.bed} 100 {params.asm}.fai | bgzip -c > {output.bins}
tabix -C {output.bins}
"""

rule GetMeanCoverage:
    input:
        bins="hmm/coverage.bins.bed.gz",
    output:
        avg="hmm/mean_cov.txt",
    shell:"""
zcat {input.bins} | awk 'BEGIN{{OFS="\\t";c=0;sum=0;}} sum=sum+$4;c=c+1;END{{print sum/c;}}' | tail -1> {output.avg}
"""

rule RunVitter:
    input:
        avg="hmm/mean_cov.txt",
        bins="hmm/coverage.bins.bed.gz",
    output:
        cov=protected("hmm/{contig}.viterout.txt"),
    params:
        sd=SD,
        contig_prefix="{contig}",
        scaler=config['scaler'],
        epsi=config['epsi']
    shell:"""
mean=$(cat {input.avg})
touch {output.cov}
tabix {input.bins} {wildcards.contig} | cut -f 4 | \
{params.sd}/viterbi  /dev/stdin $mean hmm/{params.contig_prefix} {params.scaler} {params.epsi}

"""

rule orderVitter:
    input:
        CopyNumber="hmm/{contig}.viterout.txt",
        bins="hmm/coverage.bins.bed.gz",
    output:
        cn="hmm/copy_number.{contig}.bed",
    shell:"""

paste <(tabix  {input.bins} {wildcards.contig} ) <(cat {input.CopyNumber}  ) > {output.cn}

"""

rule combineVitter:
    input:
        allCopyNumberBED=expand("hmm/copy_number.{contig}.bed", contig=contigs),
    output:
        allCN="hmm/copy_number.tsv",
    run:
        import sys
        sys.stderr.write("Writing " + str(output.allCN) + "\n")
        cn=open(output.allCN,'w')
        for fn in input.allCopyNumberBED:
            sys.stderr.write("writing " + str(fn) + "\n")
            f=open(fn)
            l=f.readlines()
            cn.write("".join(l))
        cn.close()



rule PlotBins:
    input:
        allCN="hmm/copy_number.tsv",
        #aln=config["aln"],
        avg="hmm/mean_cov.txt",
    output:
        plot="hmm/{prefix_bam}.noclip.pdf",
    params:
        sd=SD,
        genome_prefix="{prefix_bam}",
    shell:"""
touch {output.plot}
#plot every 50000 points ~ 5MB
#Rscript {params.sd}/plot.HMM.noclip.R {input.allCN} {params.genome_prefix} 50000 {input.avg}
touch {output.plot}
"""
