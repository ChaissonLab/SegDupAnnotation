import os
import tempfile
import subprocess
import os.path


# Config
configfile: "/project/mchaisso_100/projects/HPRC/sd_analysis.json"


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

assm="/project/mchaisso_100/shared/references/hg38_noalts/hg38.no_alts.fasta"




ASM=SD+"/hmcnc/HMM/annotation/hg38.fa.fai"
REP=SD+"/hmcnc/HMM/annotation/repeatMask.merged.bed"
GEN=SD+"/hmcnc/HMM/annotation/gencode.gene.bed"


fai= open(ASM)

contigs = [l.split()[0].strip().replace("|","_") for l in fai]



ep=config["epsi"]

config["ref_windows"]="/project/mchaisso_100/projects/HPRC/ref_windows.bed"

tempp=config['temp']
if "temp2" not in config:
    config["temp2"] = config["temp"]
    
if config['temp2']!="":
    tempp=config['temp2']




rule all:
    input:
        fai=assm+".fai",
        bed="hmm_ref/cov.bed",
        no_subreadbed="hmm_ref/cov.no_subread.bed",
        bins="hmm_ref/coverage.bins.bed.gz",
        avg="hmm_ref/mean_cov.txt",
        vitterout=expand("hmm_ref/{ctg}.viterout.txt", ctg=contigs),
        cn=expand("hmm_ref/copy_number.{ctg}.bed", ctg=contigs),
        allCN="hmm_ref/copy_number.tsv",
        dups="hmm_ref/collapsed_duplications.split.bed",
        don="Rhmm.done"
# Simple preprocessing, make sure there is an index on the assembly.
#




rule MakeCovBed:
    input:
        bam="ref_aligned.bam",
    output:
        bed="hmm_ref/cov.bed",
    params:
        sd=SD,
    shell:"""
samtools view -q 10 -F 2304 -@ 3 {input.bam} | {params.sd}/hmcnc/HMM/samToBed /dev/stdin/ --useH --flag   > {output.bed}
"""

if config['index_params']==" -CLR":    
    rule FilterSubreads:
        input:
            bed="hmm_ref/cov.bed",
        output:
            covbed="hmm_ref/cov.no_subread.bed",
        params:
            sd=SD,
        shell:"""
    {params.sd}/RemoveRedundantSubreads.py --input {input.bed} |sort -k1,1 -k2,2n > {output.covbed}
    """
else:
    rule FilterSubreads0:
        input:
            bed="hmm_ref/cov.bed",
        output:
            covbed="hmm_ref/cov.no_subread.bed",
        params:
            sd=SD,
        shell:"""
cd hmm_ref; ln -s cov.bed cov.no_subread.bed
    """


 
rule MakeIntersect:
    input:
        windows=config["ref_windows"],  ## add to json
        bed="hmm_ref/cov.no_subread.bed",
    output:
        bins="hmm_ref/coverage.bins.bed.gz",
    params:
        asm=assm,
        sd=SD,
    shell:"""
intersectBed -sorted -c -a {input.windows} -b {input.bed}| bgzip -c > {output.bins}
tabix -C {output.bins}
"""


rule GetMeanCoverage:
    input:
        bins="hmm_ref/coverage.bins.bed.gz",
    output:
        avg="hmm_ref/mean_cov.txt",
    shell:"""
zcat {input.bins} | awk 'BEGIN{{OFS="\\t";c=0;sum=0;}} sum=sum+$4;c=c+1;END{{print sum/c;}}' | tail -1> {output.avg}
"""

rule RunVitter:
    input:
        avg="hmm_ref/mean_cov.txt",
        bins="hmm_ref/coverage.bins.bed.gz",
    output:
        cov="hmm_ref/{contig}.viterout.txt",
    params:
        sd=SD,
        contig_prefix="{contig}",
        scaler=config['scaler'],
        epsi=config['epsi'],
    shell:"""
mean=$(cat {input.avg})
touch {output.cov}
tabix {input.bins} {wildcards.contig} | cut -f 4 | \
{params.sd}/hmcnc/HMM/viterbi  /dev/stdin $mean hmm_ref/{params.contig_prefix} {params.scaler} {params.epsi} 0

"""

rule orderVitter:
    input:
        CopyNumber="hmm_ref/{contig}.viterout.txt",
        bins="hmm_ref/coverage.bins.bed.gz",
    output:
        cn="hmm_ref/copy_number.{contig}.bed",
    shell:"""

paste <(tabix  {input.bins} {wildcards.contig} ) <(cat {input.CopyNumber}  ) > {output.cn}

"""

rule combineVitter:
    input:
        allCopyNumberBED=expand("hmm_ref/copy_number.{contig}.bed", contig=contigs),
    output:
        allCN="hmm_ref/copy_number.tsv",
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



## prefix bam
rule PlotBins:
    input:
        allCN="hmm_ref/copy_number.tsv",
        #aln=config["aln"],
        avg="hmm_ref/mean_cov.txt",
    output:
        plot="hmm_ref/ref_genome.noclip.pdf",
    params:
        sd=SD,
    shell:"""
#touch {output.plot}
#plot every 50000 points ~ 5MB
Rscript {params.sd}/hmcnc/HMM/plot.HMM.noclip.R {input.allCN} hmm_ref/ref_genome 50000 {input.avg}
#touch {output.plot}
"""


rule ConvertHMMCopyNumberToCollapsedDuplications:
    input:
        bed="hmm_ref/copy_number.tsv"
    output:
        dups="hmm_ref/collapsed_duplications.bed"
    resources:
        load=1
    shell:"""
cat {input.bed} | awk '{{ if ($5 > 2) print;}}' > {output.dups}
"""

rule MakeCoverageBins:
    input:
        cb="hmm_ref/collapsed_duplications.bed"
    output:
        cbcol="hmm_ref/collapsed_duplications.bed.collapse",
        s="hmm_ref/pre.collapsed_duplications.split.bed"
    params:
        sd=SD
    resources:
        load=1
    shell:"""
bedtools merge -i {input.cb} -c 5 -o collapse > {output.cbcol}
{params.sd}/SplitCoverageBins.py {output.cbcol} | bedtools merge -c 4,4 -o min,max > {output.s}
"""



rule repeatMask:
    input:
        call="hmm_ref/DUPcalls.copy_number.tsv",
    output:
        maskedcall="hmm_ref/DUPcalls.masked_CN.tsv",
        inter="hmm_ref/DUPcalls.rep_int.bed",
    params:
        rep=REP,
        sd=SD,
    shell:"""

    intersectBed -wa -wb -a <( awk '$4>2' {input.call} |cut -f 1-3) -b {params.rep} |sort -k1,1 -k2,2n | python {params.sd}/hmcnc/HMM/repeatMask.py | groupBy -g 1,2,3,10 -c 9| awk 'BEGIN{{OFS="\t"}} $6=$5/$4' > {output.inter}

intersectBed -v -a <( awk '$4>2' {input.call} |cut -f 1-3) -b {params.rep} |awk 'BEGIN{{OFS="\t"}}{{print$1,$2,$3,$3-$2,0,0}}'>>{output.inter}

intersectBed -wa -wb -a <(awk '$6<0.8' {output.inter} ) -b <( awk '$4>2' {input.call} ) |awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$10,$6;}}'|awk '$3-$2>15000' > {output.maskedcall}
"""



rule Postcn3:
    input:
        s="hmm_ref/pre.collapsed_duplications.split.bed"
    output:
        pre="hmm_ref/pre_cn3.txt",
        reg="hmm_ref/cn3_region.txt",
        nf="hmm_ref/cn3.nucfreq.bed.gz"
    params:
        sd=SD,
        bam="ref_aligned.bam",
        asm=assm,
        temp=config['temp']
    resources:
        load=1
    shell:"""
awk ' {{if ($4==$5 && $4==3) print ;}}' {input.s} | sort -k1,1 -k2,2n > {output.pre}

awk '{{print $1":"$2"-"$3}}' {output.pre} > {output.reg}

{params.sd}/bamToFreq {params.bam} {output.reg} {params.asm}| awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$2+1,$3,$4,$5,$6; }} ' |  bgzip -c > {output.nf}

tabix -C {output.nf}


"""


rule lrt:
    input:
        nf="hmm_ref/cn3.nucfreq.bed.gz",
        reg="hmm_ref/cn3_region.txt",
    output:
        post="hmm_ref/cn3/post_cn3.bed",
    params:
        sd=SD,
    resources:
        load=1
    shell:"""
    rm -f {output}
    echo "filtering cn3"
    for r in ` cat {input.reg} `;do
        echo $r
        tabix {input.nf} $r |  python {params.sd}/het_check.ini.py -r $r | tr ":-" "\\t" >> {output}
    done
    
"""



rule filterCN3:
    input:
        post="hmm_ref/cn3/post_cn3.bed",
        s="hmm_ref/pre.collapsed_duplications.split.bed"
    output:
        ss="hmm_ref/collapsed_duplications.split.bed"
    resources:
        load=1
    shell:"""
intersectBed -v -a {input.s} -b <( cat {input.post} |grep fail) > {output.ss} 
    """



rule callSummary:
    input:
        call="hmm_ref/DUPcalls.masked_CN.tsv",
    output:
        sumcall="hmm_ref/CallSummary.tsv",
        Vsumcall="hmm_ref/CallSummary.verbose.tsv",
    shell:"""
sort -k4,4n {input.call}|awk 'BEGIN{{OFS="\t"}} $6=$3-$2' -|groupBy -g 4 -c 4,6 -o count,mean>{output.sumcall}
sort -k1,1 -k4,4n {input.call}|awk 'BEGIN{{OFS="\t"}} $6=$3-$2' -|groupBy -g 1,4 -c 4,6 -o count,mean> {output.Vsumcall}
    """


rule GeneCount:
    input:
        call="hmm_ref/DUPcalls.masked_CN.composite.tsv",
    output:
        geneCount="hmm_ref/DUP.gene_count.bed",
    params:
        gen=GEN,
    shell:"""
        intersectBed -F 1 -wb -wa -a {input} -b {params.gen} | groupBy -g 1,2,3,6 -c 5,10,10 -o collapse,collapse,count > {output}
        """



rule Done:
    input:
        bam="ref_aligned.bam",
        s="hmm_ref/collapsed_duplications.split.bed",
    output:
        don="Rhmm.done"
    shell:"""
touch {output}
    """

