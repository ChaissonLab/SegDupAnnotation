import os
import tempfile
import subprocess
import os.path


# Config
configfile: "sd_analysis.json"

RD = config["scr"]

ASM=RD+"/annotation/hg38.fa.fai"
REP=RD+"/annotation/repeatMask.merged.bed"
GEN=RD+"/annotation/gencode.gene.bed"

#config("hmm_caller.json")


#configfile: config["configfile"]
if config["outdir"]==".":
    outdir="hmm"
else:
    outdir=config["outdir"]+"/hmm"

fai= open(ASM)

contigs = [l.split()[0].strip().replace("|","_") for l in fai]

#fofn

bamt =  config["bam"]

bamm = bamt.split("/")[-1]
prefix_bam = os.path.splitext(bamm)[0]


ep=config["epsi"]


cov=config["coverage"]
#mask=config["repeatMask"]

sub=config["subread"]
mq=config["mq"]
ml=config["minL"]


tempp=config['temp']
if "temp2" not in config:
    config["temp2"] = config["temp"]
    
if config['temp2']!="":
    tempp=config['temp2']


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

assm="/project/mchaisso_100/shared/references/hg38_noalts/hg38.no_alts.fasta"



bamFiles={f.split("/")[-1]: f for f in config["reads_bam"] }




rule all:
    input:
        fai=assm+".fai",
        bam=config["bam"],
        bed="hmm/cov.bed",
        no_subreadbed="hmm/cov.no_subread.bed",
        bins="hmm/coverage.bins.bed.gz",
        avg="hmm/mean_cov.txt",
        vitterout=expand("hmm/{ctg}.viterout.txt", ctg=contigs),
        cn=expand("hmm/copy_number.{ctg}.bed", ctg=contigs),
        allCN="hmm/copy_number.tsv",
#        plot=expand("hmm/{bm}.noclip.pdf",bm=prefix_bam),




        dups="collapsed_duplications.bed",

# Simple preprocessing, make sure there is an index on the assembly.
#


rule MakeFaiLinkOrig:
    input:
        asm=assm,
    output:
        orig="assembly.hg38.fa",
        fai=assm+".fai"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
ln -s {input.asm} ./{output.orig}
samtools faidx {output.orig}
"""


rule IndexGenome:
    input:
        ref=assm,
    output:
        gli=assm+".gli"
    params:
        sd=SD,
        grid_opts=config["grid_large"],
        index_params=config["index_params"]
    shell:"""
lra index {input.ref} {params.index_params}
"""

def GetBam(f):
    return f.split("/")[-1]

#
# Map individual bams separately
#
    
rule AlignBam:
    input:
        bam=lambda wildcards: bamFiles[wildcards.base],
        gli=assm+".gli"
    output:
        aligned="ref_aligned/{base}.bam"
    params:
        sd=SD,
        ref=assm,
        grid_opts=config["grid_large"],
        temp=config["temp"],
        mapping_params=config["mapping_params"]
    resources:
        load=16
    shell:"""

#{params.sd}/Cat.sh {input.bam} | ./home1/mchaisso/projects/LRA/lra/lra align {params.ref} - -t 16 -p s {params.mapping_params} | \
 #  samtools sort -T {params.temp}/asm.$$ -m2G -o {output.aligned}

{params.sd}/Cat.sh {input.bam} | minimap2 {params.ref} - -t 16 -a  | \
   samtools sort -T {params.temp}/asm.$$ -m2G -o {output.aligned} 

"""

rule MergeBams:
    input:
        aln=expand("ref_aligned/{b}.bam", b=bamFiles.keys())
    output:
        bam="ref_aligned.bam",
     params:
        grid_opts=config["grid_medium"]
    resources:
        load=2
    shell:"""
samtools merge {output.bam} {input.aln} -@2
"""

rule IndexBam:
    input:
        bam="ref_aligned.bam"
    output:
        bai="ref_aligned.bam.bai"
    resources:
        load=2
    params:
        grid_opts=config["grid_medium"]
    shell:"""
samtools index -@2 {input.bam}
"""


rule MakeCovBed:
    input:
        bam="ref_aligned.bam",
    output:
        bed="hmm/cov.bed",
    params:
        sd=SD,
    shell:"""
mkdir -p hmm
samtools view -q 10 -F 2304 -@ 3 {input.bam} | {params.sd}/hmcnc/HMM/samToBed /dev/stdin/ --useH --flag   > {output.bed}
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
        windows=config["ref_window"],  ## add to json
        bed="hmm/cov.no_subread.bed",
    output:
        bins="hmm/coverage.bins.bed.gz",
    params:
        asm=assm,
        sd=SD
    shell:"""
intersectBed -sorted -c -a {input.windows} -b {input.covbed}| bgzip -c > {output.bins}
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
        cov="hmm/{contig}.viterout.txt",
    params:
        sd=SD,
        contig_prefix="{contig}",
        scaler=config['scaler'],
        epsi=config['epsi']
    shell:"""
mean=$(cat {input.avg})
touch {output.cov}
tabix {input.bins} {wildcards.contig} | cut -f 4 | \
{params.sd}/hmcnc/HMM/viterbi  /dev/stdin $mean hmm/{params.contig_prefix} {params.scaler} {params.epsi} 0

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



## prefix bam
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
#touch {output.plot}
#plot every 50000 points ~ 5MB
Rscript {params.sd}/hmcnc/HMM/plot.HMM.noclip.R {input.allCN} hmm/{params.genome_prefix} 50000 {input.avg}
#touch {output.plot}
"""


rule ConvertHMMCopyNumberToCollapsedDuplications:
    input:
        bed="hmm/copy_number.tsv"
    output:
        dups="hmm/collapsed_duplications.bed"
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
cat {input.bed} | awk '{{ if ($5 > 2) print;}}' > {output.dups}
"""

rule MakeCoverageBins:
    input:
        cb="hmm/collapsed_duplications.bed"
    output:
        cbcol="hmm/collapsed_duplications.bed.collapse",
        s="hmm/pre.collapsed_duplications.split.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
bedtools merge -i {input.cb} -c 5 -o collapse > {output.cbcol}
{params.sd}/SplitCoverageBins.py {output.cbcol} | bedtools merge -c 4,4 -o min,max > {output.s}
"""



rule repeatMask:
    input:
        call="{outdir}/hmm_ref/DUPcalls.copy_number.tsv",
    output:
        maskedcall="{outdir}/hmm_ref/DUPcalls.masked_CN.tsv",
        inter="{outdir}/hmm_ref/DUPcalls.rep_int.bed",
    params:
        rep=REP,#config['repeatMask'],
        mnl=config["minL"],
        rd=RD,
    shell:"""

    intersectBed -wa -wb -a <( awk '$4>2' {input.call} |cut -f 1-3) -b {params.rep} |sort -k1,1 -k2,2n | python {params.rd}/repeatMask.py | groupBy -g 1,2,3,10 -c 9| awk 'BEGIN{{OFS="\t"}} $6=$5/$4' > {output.inter}

intersectBed -v -a <( awk '$4>2' {input.call} |cut -f 1-3) -b {params.rep} |awk 'BEGIN{{OFS="\t"}}{{print$1,$2,$3,$3-$2,0,0}}'>>{output.inter}

intersectBed -wa -wb -a <(awk '$6<0.8' {output.inter} ) -b <( awk '$4>2' {input.call} ) |awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$10,$6;}}'|awk '$3-$2>{params.mnl}' > {output.maskedcall}
"""



rule Postcn3:
    input:
        s="hmm/pre.collapsed_duplications.split.bed"
    output:
        pre="hmm/pre_cn3.txt",
        reg="hmm/cn3_region.txt",
        nf="hmm/cn3.nucfreq.bed.gz"
    params:
        grid_opts=config["grid_blat"],
        sd=SD,
        bam=config['bam'],
        asm=assembly,
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
        nf="hmm/cn3.nucfreq.bed.gz",
        reg="hmm/cn3_region.txt",
     #   done="getPos.done",
       # ps=lambda wildcards: pos[wildcards.p],
    output:
        post="hmm/cn3/post_cn3.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SD,
        #ps="{p}",
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
        post="hmm/cn3/post_cn3.bed",
        #expand("cn3/post_cn3.{p}.bed",p=lambda wildcards: getPos("cn3_region.txt")),
        s="hmm/pre.collapsed_duplications.split.bed"
    output:
        ss="collapsed_duplications.split.bed"
    resources:
        load=1
    params:
        grid_opts=config["grid_small"],
    shell:"""
intersectBed -v -a {input.s} -b <( cat {input.post} |grep fail) > {output.ss} 
    """



rule callSummary:
    input:
        call="{outdir}/hmm_ref/DUPcalls.masked_CN.tsv",
    output:
        sumcall="{outdir}/hmm_ref/CallSummary.{ep}.tsv",
        Vsumcall="{outdir}/hmm_ref/CallSummary.verbose.{ep}.tsv",
    shell:"""
sort -k4,4n {input.call}|awk 'BEGIN{{OFS="\t"}} $6=$3-$2' -|groupBy -g 4 -c 4,6 -o count,mean>{output.sumcall}
sort -k1,1 -k4,4n {input.call}|awk 'BEGIN{{OFS="\t"}} $6=$3-$2' -|groupBy -g 1,4 -c 4,6 -o count,mean> {output.Vsumcall}
    """


rule GeneCount:
    input:
        call="{outdir}/hmm_ref/DUPcalls.masked_CN.composite.tsv",
    output:
        geneCount="{outdir}/hmm_ref/DUP.gene_count.bed",
    params:
        gen=GEN,
    shell:"""
        intersectBed -F 1 -wb -wa -a {input} -b {params.gen} | groupBy -g 1,2,3,6 -c 5,10,10 -o collapse,collapse,count > {output}
        """




