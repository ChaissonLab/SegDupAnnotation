import os
import tempfile
import subprocess
import os.path
# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "sd_analysis.json"

assembly=config["asm"]
geneModel = config["genemodel"].keys()

RNADatasets = config["RNAseq"].keys()
spliced=[ "multi", "single"]

if "TMPDIR" not in os.environ:
    print("ERROR, TMPDIR must be set")
    sys.exit(1)

rule all:
    input:
        fai=assembly+".fai",
        wm_db="wmdb",
        wm_intv="wm_mask_intervals",
        masked="assembly.masked.fasta",
        sedef="sedef_out/final.bed",
        sedef_sorted="sedef_out/final.sorted.bed",
        sdcount="sedef_out2/final.sorted.bed.filt.count",
        sd2count="sedef_out2/final.low_copy.bed",
        sedef_filtered="sedef_out/final.sorted.filtered.bed",
        pairNotMasked="sedef_out2/final.low_copy.not_masked_pair.bed",
        pairMasked="sedef_out2/final.low_copy.masked_pair.bed",
        genesInResolvedDups="genes_in_resolved_dups.bed",
        oneIsoform="genes_in_resolved_dups.one_isoform.bed",
        realignedOneIsoform="genes_in_resolved_dups.one_isoform.bed.bam",
        b12="genes_in_resolved_dups.one_isoform.bed.bam.bed12",
        splitAndSpliced=expand("genes_in_resolved_dups.one_isoform.{sp}.bed", sp=spliced),
        alignedIsoforms=expand("identity.{sp}.bed", sp=spliced),
        missed_repeats="sedef_out/missed.bed",
        sedef_filt="sedef_out/final.sorted.bed.filt",
        qf="sedef_out/final.sorted.bed.qual",
        repmasked="assembly.repeat_masked.fasta",
        dups="collapsed_duplications.bed",
        hmmCopyNumber="hmm/copy_number.tsv",
        rnabam=expand("{data}.mapped.bam",data=geneModel),
        rnabambed=expand("{data}.mapped.bam.bed12",data=geneModel),
        rnabambed6=expand("{data}.mapped.bam.bed6",data=geneModel),
        sup=expand("{data}.mapped.bam.bed6.rnaseq-sup",data=geneModel),
        supdup=expand("{data}.mapped.bam.bed12.dups.sup",data=geneModel),
        colrnabed=expand("{data}.mapped.bam.bed12.dups",data=geneModel),
        resdup=expand("{data}.mapped.resolved_dups.bed",data=geneModel),
        gencodeRes="gencode.mapped.multicopy.bed",
        gencodeResSup="gencode.mapped.multicopy.bed.supported",
        duplicationSummary="gencode.combined-duplicated-genes.tsv",
        summaryGencodeResSup="gencode.resolved-duplications.tsv",
        gencodeResIDup="gencode.mapped.multicopy.in_duplication.bed",
        resdupsup=expand("{data}.mapped.resolved_dups.bed.sup",data=geneModel),
        combineMasked="assembly.union_masked.fasta",
        RNAseq=expand("RNAseq/{dataset}.bam", dataset=list(config["RNAseq"].keys())),
        RNAseqCov=expand("RNAseq/{dataset}.bam.cov", dataset=list(config["RNAseq"].keys())),
        combinedCov="RNAseq/combined.bed",
        IsoSeq=expand("IsoSeq/{dataset}.bam", dataset=list(config["IsoSeq"].keys())),
        filt="sedef_out/final.sorted.bed.qualfilt",
        counted="sedef_out/counted.bed",
        countedhc="sedef_out/counted.bed.high_copy",
        rep="sedef_out/counted.bed.rep",
        rgn="sedef_out/counted.bed.rep.rgn",
        ident="sedef_out/counted.bed.rep.rgn.ident",
        countedall="sedef_out/final.sorted.counted.bed",
        low_copy="sedef_out/final.low_copy.sorted.bed",
        low_copy_cov="depth_over_dup.tsv",
        low_copy_cov_filt="depth_over_dup.filt.tsv",
        average_cov="average_coverage.txt",
        sdMask="assembly.repeat_masked.sd.fasta",
        sedef2="sedef_out2/final.bed",
        sedef2Sorted="sedef_out2/final.sorted.bed",
        sedef2SortedFilt="sedef_out2/final.sorted.bed.filt",
        asmMask=expand("{asm}.count_masked", asm=["assembly.orig.fasta", "assembly.masked.fasta", "assembly.repeat_masked.fasta", "assembly.repeat_masked.sd.fasta"]),
        uniqueDupGenes="gencode.mapped.bam.bed12.dups.unique",
        uniqueDupGenesCN="gencode.mapped.bam.bed12.dups.unique.cn"

rule GetAverageCoverage:
    input:
        bam=config["bam"]
    output:
        avg="average_coverage.txt",
    params:
        sd=SD
    shell:"""
samtools depth {input.bam} | awk '{{ if (NR%1000==0)  print;}}' | head -3000 | cut -f 3  | {params.sd}/stats > average_coverage.txt
"""

rule WriteCoverageOverDups:
    input:
        dups="sedef_out2/final.low_copy.not_masked_pair.bed",
        bam=config["bam"]
    params:
        sd=SD,
    output:
        dod="depth_over_dup.tsv"
    shell:"""
{params.sd}/MeasureDepthOfResolvedDuplications.sh {input.dups} {input.bam} > {output.dod}
"""

rule FilterDepthOverDups:
    input:
        dod="depth_over_dup.tsv",
        avg="average_coverage.txt"
    output:
        dodf="depth_over_dup.filt.tsv"
    shell:"""
avg=`cut -f 1 {input.avg}`
cat {input.dod} | awk -v avg=$avg '{{ if (($3 < avg*0.3) || ($4 < avg*0.3)) {{  print;}} }}' > {output.dodf}
"""

rule RemapBed:
    input:
        oneIsoform="genes_in_resolved_dups.one_isoform.bed",
        asm="assembly.orig.fasta",
        refGenes=config["genemodel"]["gencode"],
        refGeneBam="gencode.mapped.bam"
    output:
        realignedOneIsoform="genes_in_resolved_dups.one_isoform.bed.bam",
    params:
        sd=SD
    shell:"""
{params.sd}/RemapRegion.sh {input.oneIsoform} {input.asm} {input.refGenes} {input.refGeneBam} {output.realignedOneIsoform}
"""
   
rule RemappedBedToSummaryTable:
    input:
        realignedOneIsoform="genes_in_resolved_dups.one_isoform.bed.bam",
    output:
        b12="genes_in_resolved_dups.one_isoform.bed.bam.bed12",
        names="genes_in_resolved_dups.one_isoform.bed.bam.name.bed12",
        counts="genes_in_resolved_dups.one_isoform.bed.bam.name.counts",
    shell:"""
bedtools bamtobed -bed12 -i {input.realignedOneIsoform} > {output.b12}
bedtools bamtobed -bed12 -i {input.realignedOneIsoform} | awk '{{ split($4,n,"|"); $4=n[6]; print $0; }}' | tr " " "\t" > {output.names}
bedtools groupby -g 4 -c 4 -i {output.names} -o count > {output.counts}
"""


rule GetUniqueGencodeUnresolvedDupGenesCN:
    input:
        uniqueDupGenes="gencode.mapped.bam.bed12.dups.unique"
    output:
        uniqueDupGenesCN="gencode.mapped.bam.bed12.dups.unique.cn"
    params:
        grid_opts=config["grid_small"],
    shell:"""
cut -f 4,16 {input.uniqueDupGenes} > {output.uniqueDupGenesCN}
"""
    
rule GetUniqueGencodeUnresolvedDupGenes:
    input:
        bed="gencode.mapped.bam.bed12.dups"
    output:
        unique="gencode.mapped.bam.bed12.dups.unique"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
cat {input.bed} |  {params.sd}/WriteNonOverlapping.py  > {output.unique}
"""
   

rule GetGencodeMulticopy:
    input:
        bed="gencode.mapped.bam.bed12",
    output:
        gencodeRes="gencode.mapped.multicopy.bed"
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:"""
cut -f 4 {input.bed} | sort | uniq -c | awk '{{ if ($1 > 1) print $2;}}' > {input.bed}.multimap
{params.sd}/FindDuplicatedResolved.py {input.bed} {input.bed}.multimap  | \
 awk '{{ print $3-$2"\\t"$0; }}' | \
 sort -k5,5 -k1,1nr | \
 cut -f 2- | \
 {params.sd}/FilterSimilarLengths.py > {output.gencodeRes} 
"""

rule GetGencodeMappedInDup:
    input:
        gencodeRes="gencode.mapped.multicopy.bed",
        dups="sedef_out2/final.sorted.bed"
    output:
        inDup="gencode.mapped.multicopy.in_duplication.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.dups} | awk '{{ if ($8 < 10) print;}}' > {input.dups}.ten
bedtools intersect -a {input.gencodeRes} -b {input.dups}.ten -u -f 1 > {output.inDup}
"""
   

rule GetSupportedMulticopy:
    input:
        bed="gencode.mapped.multicopy.bed",
        sup="gencode.mapped.bam.bed6.rnaseq-sup"
    output:
        sup="gencode.mapped.multicopy.bed.supported",
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -u -a {input.bed}  -b {input.sup} > {output.sup}
"""


rule SummarizeMulticopy:
    input:
        sup="genes_in_resolved_dups.one_isoform.bed"
    output:
        summary="gencode.resolved-duplications.tsv",
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.sup} | \
  awk '{{ print $13"\\t"$14"\\t"$15"\\t"$4; print $16"\\t"$17"\\t"$18"\\t"$4;}}' | \
  sort -k4,4 -k1,1 -k2,2n  | \
  bedtools groupby -g 1-4 -o first -full -c 4 | cut -f 4 | uniq -c | awk '{{ print $2"\\t"$1"\\tresolved"; }}' > {output.summary}
"""

rule UnionMasked:
    input:
        orig="assembly.orig.fasta",
        wm="assembly.masked.fasta",
        rm="assembly.repeat_masked.fasta"
    output:
        comb="assembly.union_masked.fasta"
    params:
        sd=SD,
        grid_opts=config["grid_medium"]
    shell:"""
{params.sd}/comask {output.comb} {input.orig} {input.wm} {input.rm}
"""

rule RunDepthHmm:
    input:
        bam=config["bam"],
        asm="assembly.orig.fasta"
    output:
        vo="hmm/copy_number.tsv"
    params:
        grid_opts=config["grid_large"],
        sd=SD
    shell:"""
snakemake --nolock -p -s {params.sd}/hmm_caller.vert.snakefile -j 20
"""

rule SelectDups:
    input:
        bed="hmm/copy_number.tsv"
    output:
        dups="collapsed_duplications.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.bed} | awk '{{ if ($5 > 3) print;}}' > {output.dups}
"""


#
# Find all the full-length genes that overlap a resolved duplication.
#
rule FindResolvedDupliatedGenes:
    input:
        rnabed="{data}.mapped.bam.bed12",
        sedef="sedef_out2/final.sorted.bed.filt",
    output:
        resdup="{data}.mapped.resolved_dups.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
bedtools intersect -a {input.rnabed} -b {input.sedef} -loj -f 1 | awk '{{ if ($13 != ".") print; }}' | \
{params.sd}/RefSeqToName.py | bedtools sort >  {output.resdup}
"""



rule Bed12ToBed6:
    input:
        pre="{base}.bed12"
    output:
        post="{base}.bed6"
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools bed12tobed6 -i {input.pre} > {output.post}
"""

rule GetRNASeqSupportedTranscript:
    input:
        bed6="{base}.bed6",
        comb="RNAseq/combined.bed"
    output:
        sup="{base}.bed6.rnaseq-sup"
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -a {input.bed6} -b {input.comb} -wao | bedtools sort > {output.sup}
"""

rule CombineGenesWithCollapsedDups:
    input:
        rnabed="{data}.mapped.bam.bed12",
        dups="collapsed_duplications.split.bed",
        asm="assembly.union_masked.fasta"
    output:
        rnabedout="{data}.mapped.bam.bed12.dups",
    params:
        grid_opts=config["grid_small"],
        sd=SD, 
    shell:"""
samtools faidx {input.asm}
bedtools intersect -f 1 -g {input.asm}.fai -sorted -loj -a {input.rnabed} -b {input.dups} | awk '{{ if ($13 != ".") print;}}'| \
{params.sd}/RefSeqToName.py | \
bedtools sort >  {output.rnabedout}
"""

rule SummarizeGenesWithCollapsedDups:
    input:
        genes="gencode.mapped.bam.bed12.dups",
        sup="gencode.mapped.bam.bed12.dups.sup"
    output:
        coll="gencode.collapsed.tsv"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
cat {input.sup} | \
 {params.sd}/FilterIsoformsFromCollapse.py | \
cut -f 4,16 |  awk '{{ print $1"\\t"$2"\\tcollapse";}}' > {output.coll}
"""

rule CombineDuplicatedGenes:
    input:
        coll = "gencode.collapsed.tsv",
        res  = "gencode.resolved-duplications.tsv",
    output:
        summary = "gencode.combined-duplicated-genes.tsv"
    params:
        grid_opts = config["grid_small"]
    shell:"""
cat {input.coll} {input.res} | sort > {output.summary}
"""

rule GetSupportedDups:
    input:
        dup="{data}.mapped.bam.bed12.dups",
        rna="{data}.mapped.bam.bed6.rnaseq-sup"
    output:
        supdup="{data}.mapped.bam.bed12.dups.sup",
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -sorted -wao -a {input.dup} -b {input.rna} | bedtools groupby -c 1 -o first -full > {output.supdup}
"""

rule GetSupportedResolvedDups:
    input:
        dup="{data}.mapped.resolved_dups.bed",
        rna="{data}.mapped.bam.bed6.rnaseq-sup"
    output:
        supdup="{data}.mapped.resolved_dups.bed.sup",
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -sorted -wao -a {input.dup} -b {input.rna} | bedtools groupby -c 1 -o first -full > {output.supdup}
"""


#rule CountDuplicatedGenes:
#    input:
#        rnabed=
#
rule LinkOrig:
    input:
        asm=config["asm"]
    output:
        orig="assembly.orig.fasta"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
ln -s {input.asm} ./{output.orig}
"""
    
rule CountMaskedAsm:
    input:
        asm="{assembly}"
    output:
        count="{assembly}.count_masked"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
cat {input.asm} | {params.sd}/nl > {output.count}
"""
        
rule CountedToHC:
    input:
        counted="sedef_out/counted.bed",
    output:
        countedhc="sedef_out/counted.bed.high_copy",
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.counted} | awk '{{ if ($4 > 20) print;}}' > {output.countedhc}
"""

rule MaskBySedef:
    input:
        asm="assembly.union_masked.fasta", 
        bed="sedef_out/counted.bed.high_copy"
    output:
        rmsk="assembly.repeat_masked.sd.fasta"
    params:
        grid_opts=config["grid_large"],
        sd=SD
    shell:"""
{params.sd}/bemask {input.asm} {input.bed} {output.rmsk}
samtools faidx {output.rmsk}
"""

rule RerunSedef:
    input:
        asm="assembly.repeat_masked.sd.fasta"
    output:
        sd2="sedef_out2/final.bed"
    params:
        grid_opts=config["grid_large"] + " --nodelist=\"hpc[4574-4577]\""
    shell:"""
export PATH=$PATH:$mchaisso/software/sedef/
sedef.sh {input.asm} -f -j 16 -t {input.asm}.translate.fa -o sedef_out2
"""

rule CountRepeatedEntriesAll:
    input:
        missed="sedef_out/final.sorted.bed"
    output:
        counted="sedef_out/final.sorted.counted.bed"
    params:
        grid_opts=config["grid_medium"]
    shell:"""
bedtools intersect -loj -a {input.missed} -b {input.missed} -sorted -f 0.8 -r | cut -f 1-3 | uniq -c | awk '{{ print $2"\\t"$3"\\t"$4"\\t"$1;}}' > {output.counted}
"""

rule SortRerun:
    input:
        sd2="sedef_out2/final.bed"
    output:
        sd2s="sedef_out2/final.sorted.bed"
    params:
        sd=SD,
        grid_opts=config["grid_medium"]
    shell:"""
cat {input.sd2} | {params.sd}/AddHeaderToSedef.py | bedtools sort -header > {output.sd2s}
"""

rule FilterFinalSedef:
    input:
        sd2s="sedef_out2/final.sorted.bed"
    output:
        sd2sf="sedef_out2/final.sorted.bed.filt"
    params:
        grid_opts=config["grid_small"]
    shell:"""
bioawk -c hdr '{{ if (NR==1 || $score < 20 && $3-$2 > 5000) print; }}' < {input.sd2s} > {output.sd2sf}
"""

rule CountMatches:
    input:
        sd2sf="sedef_out2/final.sorted.bed.filt",
        unfilt="sedef_out2/final.sorted.bed"
    output:
        sdcount="sedef_out2/final.sorted.bed.filt.count"
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -loj -f 0.9 -r -sorted -a {input.sd2sf} -b {input.unfilt} | bedtools groupby -g 1-3 -c 2 -o count | awk '{{ if ($4 < 20) print;}}' > {output.sdcount}
"""

rule RetainLowCopy:
    input:
        sd2sf="sedef_out2/final.sorted.bed.filt",
        sdcount="sedef_out2/final.sorted.bed"
    output:
        sd2count="sedef_out2/final.low_copy.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -a {input.sd2sf} -b {input.sdcount} -f 1 -r -u > {output.sd2count}
"""

rule RemoveMasked:
    input:
        sd2count="sedef_out2/final.low_copy.bed",
        asm="assembly.repeat_masked.fasta"
    output:
        notmasked="sedef_out2/final.low_copy.not_masked_pair.bed",
        masked="sedef_out2/final.low_copy.masked_pair.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/FilterTooMasked.py --genome {input.asm} --sedef {input.sd2count} --pct 0.8 --maskedBed {output.masked} > {output.notmasked}.full
bedtools intersect -f 0.9 -r -v -a {output.notmasked}.full -b {output.masked} | bedtools sort | \
{params.sd}/WriteNonOverlappingPairs.py | uniq > {output.notmasked}
"""

rule FindGenesInResolvedDups:
    input:
        notmasked="sedef_out2/final.low_copy.not_masked_pair.bed",
        gencode="gencode.mapped.bam.bed12"
    output:
        gencodeInDups="genes_in_resolved_dups.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat '{input.notmasked}' | awk '{{ print $0"\t"NR; a=$1;b=$2;c=$3;d=$4;e=$5;f=$6; $1=d;$2=e;$3=f; $4=a;$5=b;$6=c; print $0"\t"NR;}}' | \
tr " " "\t" | \
bedtools intersect -loj -a {input.gencode} -b stdin -f 1  | awk '{{ if ($13 != ".") print ;}}'  > {output.gencodeInDups}
"""

rule SelectOneIsoform:
    input:
        anyIsoform="genes_in_resolved_dups.bed",
    output:
        oneIsoform="genes_in_resolved_dups.one_isoform.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/FilterMembersFromSameIsoformSet.py {input.anyIsoform}  > {output.oneIsoform}
"""
rule SplitSplicedAndSingleExon:
    input:
        oneIsoform="genes_in_resolved_dups.one_isoform.bed",
    output:
        multiAndOne=expand("genes_in_resolved_dups.one_isoform.{sp}.bed", sp=spliced),
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/SplitSplicedGenes.py {input.oneIsoform} {output.multiAndOne}
"""
    
rule AlignIsoforms:
    input:
        asm="assembly.repeat_masked.fasta",
        refGenes=config["genemodel"]["gencode"],
        dups="genes_in_resolved_dups.one_isoform.{sp}.bed"
    output:
        alignedIsoforms="identity.{sp}.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/AlignIsoforms.sh {input.dups} {input.refGenes} {input.asm} > {output.alignedIsoforms}
"""

 

rule CountRepeatedEntries:
    input:
        missed="sedef_out/missed.bed"
    output:
        counted="sedef_out/counted.bed"
    params:
        grid_opts=config["grid_medium"]
    shell:"""
bedtools intersect -loj -a {input.missed} -b {input.missed} -sorted -f 0.8 -r | cut -f 1-3 | uniq -c | awk '{{ print $2"\\t"$3"\\t"$4"\\t"$1;}}' > {output.counted}
"""

rule ExtractRepeatsFromCounted:
    input:
        counted="sedef_out/counted.bed"
    output:
        reps="sedef_out/counted.bed.rep"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.counted} | awk '{{ if ($4 > 30) print; }}' | bedtools merge -c 1 -o count > {output.reps}
"""

rule MakeRepeatRgn:
    input:
        reps="sedef_out/counted.bed.rep"
    output:
        rgn="sedef_out/counted.bed.rep.rgn"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.reps} | awk '{{ print $1":"$2"-"$3;}}' > {output.rgn}
"""

rule CountRepeatIdentity:
    input:
        reps="sedef_out/counted.bed.rep.rgn",
        asm="assembly.repeat_masked.fasta"
    output:
        ident="sedef_out/counted.bed.rep.rgn.ident"
    params:
        grid_opts=config["grid_small"],
        sd=SD,    
    shell:"""
samtools faidx {input.asm} -r {input.reps} | {params.sd}/nl > {output.ident}
"""
    
rule MapRNASeq:
    input:
        rna=lambda wildcards : config["RNAseq"][wildcards.dataset],
        asm="assembly.orig.fasta"
    output:
        aln="RNAseq/{dataset}.bam"
    params:
        grid_opts=config["grid_large"]
    shell:"""
mkdir -p RNAseq
minimap2 -t 16 -ax splice -uf -C5 {input.asm} {input.rna} | samtools view -uS - | samtools sort -T $TMPDIR/{wildcards.dataset}.$$ -m2G -o {output.aln}
samtools index {output.aln}
"""

rule RNASeqCov:
    input:
        aln="RNAseq/{dataset}.bam",
        asm="assembly.orig.fasta"
    output: 
        cov="RNAseq/{dataset}.bam.cov"
    params:
        grid_opts=config["grid_large"],
        sd=SD
    shell:"""
cat {input.asm}.fai | awk '{{ if ($2 > 20000) {{ print $1":1-"$2;}} }}'> {input.asm}.rgn

cat {input.asm}.rgn | xargs -P16 -I {{}} {params.sd}/BamToFreq.sh {input.aln} {{}} {input.asm} {output.cov}
cat {output.cov}.*.rgn | bedtools sort > {output.cov} 

"""    

rule CombineRNASeqCov:
    input:
        rna=expand("RNAseq/{dataset}.bam.cov", dataset=config["RNAseq"].keys())
    output:
        combinedCov="RNAseq/combined.bed",
    params:
        grid_opts=config["grid_medium"]
    shell:"""
cat {input.rna} | bedtools sort | bedtools merge -c 4 -o max > {output.combinedCov}
"""

rule MapIsoSeq:
    input:
        rna=lambda wildcards : config["IsoSeq"][wildcards.dataset],
        assembly="assembly.orig.fasta"
    output:
        aln="IsoSeq/{dataset}.bam" 
    params:
        grid_opts=config["grid_large"]
    shell:"""
mkdir -p IsoSeq
minimap2 -t 16 -ax splice -uf -C5 {input.assembly} {input.rna} | samtools view -uS - | samtools sort -T $TMPDIR -m2G -o {output.aln}
samtools index {output.aln}
"""

rule SimpleFilter:
    input:
        sedef="sedef_out/final.sorted.bed"
    output:
        filt="sedef_out/final.sorted.bed.qualfilt"
    params:
        grid_opts=config["grid_small"]
    shell:"""
#
# 
cat {input.sedef} | awk '{{ if ($8 < 40) print;}}' > {output.filt}

"""
    

rule FindMissedRepeats:
    input:
        sedef="sedef_out/final.sorted.bed"
    output:
        missed="sedef_out/missed.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
#
# 
cat {input.sedef} | awk '{{ if ($8 < 40) {{ a=$3-$2; b=$6-$5; if (a > b) {{ d=a-b;}} else {{ d=b-a; }} if (d < 100) print;}} }}' >{output.missed}
"""

rule FindMissedRepeatTargets:
    input:
        sedef="sedef_out/missed.bed",
        final="sedef_out/final.sorted.bed",
    output:
        target="sedef_out/missed.target.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -a {input.sedef}  -b {input.final} -loj -f 0.8 -r -sorted | awk '{{ a=$7-$6; b=$10-$9; if ($5 != "." && ( (a > b && a-b < 200) || (b > a && b-a < 200 ) ) ) {{ print $5"\\t"$6"\\t"$7; print $8"\\t"$9"\\t"$10; }}  }}' | bedtools sort | bedtools merge -c 1 -o count > {output.target}
"""

rule FilterMissedRepeats:
    input:
        sedef="sedef_out/final.sorted.bed",
        missed="sedef_out/missed.bed"
    output:
        filt="sedef_out/final.sorted.bed.filt"
    params:
        grid_opts=config["grid_blat"]
    shell:"""
bedtools intersect -v -a {input.sedef} -b {input.missed} -f 0.8 > {output.filt}
"""

rule QualFilterSedef:
    input:
        filt="sedef_out/final.sorted.bed.filt"
    output:
        qf="sedef_out/final.sorted.bed.qual"
    params:
        grid_opts=config["grid_blat"],
        sd=SD
    shell:"""
cat {input.filt} | awk '{{ if ($8 < 20) print;}}' > {output.qf}
"""
    
rule FilterPSL:
    input:
        psl="{base}.psl"
    output:
        filt="{base}.psl.filt"
    params:
        grid_opts=config["grid_blat"],
        sd=SD
    shell:"""
{params.sd}/FilterPSL.py {input.psl} 100 > {output.filt}
"""


rule SplitRNASeq:
    input:
        fa=lambda wildcards: config["genemodel"][wildcards.data]
    output:
        splitRNA=dynamic("split-rna/{data}-{id}.fa")
    params:
        grid_opts=config["grid_blat"]
    shell:"""
mkdir -p split-rna
L=`wc -l {input.fa}.fai | tr " " "\\t" | cut -f 1`
W=$(($L/100))
for i in {{0..100}}; do
  s=$(($L*($i+1)));
  c=`head -$s {input.fa}.fai | tail -$W | cut -f 1 | tr "\\n" " "`
  samtools faidx {input.fa} $c > split-rna/{wildcards.data}-$i.fa
done
"""

rule MinimapGeneModel:
    input:
        fa=lambda wildcards: config["genemodel"][wildcards.data],
        asm="assembly.orig.fasta"
    output:
        bam="{data}.minimapped.bam"
    params:
        grid_opts=config["grid_large"],
        sd=SD
    shell:"""
minimap2 -x splice -a -t 8 {input.asm} {input.fa}  | samtools sort -T $TMPDIR/tmp.$$ -o {output.bam}
samtools index {output.bam}
"""

rule FilterGeneModel:
    input:
        mmbam="{data}.minimapped.bam"
    output:
        bam="{data}.mapped.bam"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/FilterMappedLength.py {input.mmbam} any | samtools view -b -o {output.bam}
samtools index {output.bam}
"""

rule MinimapGeneModelBed:
    input:
        bam="{data}.mapped.bam"
    output:
        bed="{data}.mapped.bam.bed12"
    params:
        grid_opts=config["grid_large"]
    shell:"""
bedtools bamtobed -bed12 -i {input.bam} > {output.bed}
"""
    
#rule AlignRNASeq:
#    input:
#        fa="split-rna/{data}-{id}.fa",
#        asm="assembly.repeat_masked.fasta"
#    output:
#        psl="split-rna/{data}-{id}.fa.psl",
#    params:
#        grid_opts=config["grid_medium"]
#    shell:"""
#blat -q=rna {input.asm} {input.fa} {output.psl} 
#"""
#
rule PSLtoBed:
    input:
        aln="rna.{data}.psl"
    output:
        bed="rna.{data}.psl.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
psl2bed < {input.aln} > {output.bed}
"""

rule CombineRNASeq:
    input:
        psl=dynamic("split-rna/{data}-{id}.fa.psl")
    output:
        aln="rna.{data}.psl"
    params:
        grid_opts=config["grid_small"]
    shell:"""
head -5 split-rna/{wildcards.data}-1.fa.psl > {output.aln}
for p in {input.psl}; do
tail -n +6 $p 
done  >> {output.aln}
"""

rule FilterSedef:
    input:
        bed="sedef_out/final.sorted.bed"
    output:
        filtered="sedef_out/final.sorted.filtered.bed",
    params:
        grid_opts=config["grid_small"]
    shell:"""
bioawk -c hdr '{{ aLen=$3-$2; bLen=$6-$5; if ($score < 20) {{ print;}} }} ' > {output.filtered} < {input.bed}
"""

rule FilterHighCopyRepeats:
    input:
        comb="sedef_out/counted.bed",
        final="sedef_out/final.sorted.bed"
    output:
        low_copy="sedef_out/final.low_copy.sorted.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.comb} | awk '{{ if ($4 > 10) print }}' > {input.comb}.high_copy;
bedtools intersect -header -v -a {input.final} -b {input.comb}.high_copy -f 0.5 -r -sorted  | \
  bioawk -c hdr '{{ if (NR==1 || (($uppercaseMatches / ($end1-$start1) > 0.20) && ($uppercaseMatches / ($end2-$start2) > 0.2)))  print; }}' >  {output.low_copy}
"""

rule SortSedef:
    input:
        bed="sedef_out/final.bed"
    output:
        s="sedef_out/final.sorted.bed"
    params:
        grid_opts=config["grid_medium"]
    shell:"""
sort -k1,1 -k2,2n {input.bed} > {output.s}
"""

rule RunSedef:
    input:
        asm="assembly.union_masked.fasta"
    output:
        done="sedef_out/final.bed"
    params:
        grid_opts=config["grid_large"]  + " --nodelist=\"hpc[4574-4577]\""
    shell:"""
export PATH=$PATH:$mchaisso/software/sedef/
sedef.sh {input.asm} -j 16 -t {input.asm}.translate.fa
"""



rule SplitGenome:
    input:
        asm="assembly.masked.fasta"
    output:
        split=dynamic("split/to_mask.{region}.fasta")
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
mkdir -p split
{params.sd}/DivideFasta.py {input.asm} 10000000 split/to_mask
"""


rule MaskFasta:
    input:
        part="split/to_mask.{region}.fasta"
    output:
        mask="split/to_mask.{region}.fasta.masked"
    params:
        grid_opts=config["grid_medium"]
    shell:"""

cp {input.part} $TMPDIR/to_mask.{wildcards.region}.fasta && \
pushd $TMPDIR &&  \
RepeatMasker  -pa 4 -species whale -s -xsmall to_mask.{wildcards.region}.fasta && \
popd && \
cp $TMPDIR/to_mask.{wildcards.region}.fasta.* split/ || true
if [ ! -e split/to_mask.{wildcards.region}.fasta.masked ]; then
  cp split/to_mask.{wildcards.region}.fasta split/to_mask.{wildcards.region}.fasta.masked
fi
"""

        
rule MakeToCombine:
    input:
        mask=dynamic("split/to_mask.{region}.fasta.masked"),
    output:
        tc="to_combine.txt"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    run:
        f=open(output.tc,'w')
        f.write("\n".join(input.mask) + "\n")
        f.close()

    
rule CombineMasked:
    input:
        mask="to_combine.txt",
        asm=config["asm"]
    output: 
        masked="assembly.repeat_masked.fasta"
    params:
        grid_opts=config["grid_medium"],
        sd=SD
    shell:"""
{params.sd}/CombineMasked.py {input.mask} {input.asm}.fai {output.masked}
"""


rule MakeCoverageBins:
    input:
        cb="collapsed_duplications.bed"
    output:
        s="collapsed_duplications.split.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
bedtools merge -i {input.cb} -c 5 -o collapse > {input.cb}.collapse
{params.sd}/SplitCoverageBins.py {input.cb}.collapse | bedtools merge -c 4,4 -o min,max > {output.s}
"""

rule MakeFAI:
    input:
        asm=assembly
    output:
        fai=assembly+".fai"
    params:
       grid_opts=config["grid_small"]
    shell:"""
samtools faidx {input.asm}
"""

rule MakeWMDB:
    input:
        asm=assembly
    output:
        wm_db="wmdb"
    params:
        grid_opts=config["grid_medium"]
    shell:"""
windowmasker -mk_counts -in {input.asm} -out {output.wm_db}  || true
"""

rule MakeWMIntv:
    input:
        wm_db="wmdb",
        asm=assembly
    output:
        intv="wm_mask_intervals"
    params:
        grid_opts=config["grid_large"]
    shell:"""
windowmasker -ustat {input.wm_db} -in {input.asm} -out {output.intv}  || true
"""

rule MakeWMBed:
    input:
        intv="wm_mask_intervals"
    output:
        bed="wm_mask_intervals.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.intv} | awk '{{ if (substr($1,0,1) == ">") {{ name=substr($1,1); }} else {{ if ($3-$1 > 100) print name"\\t"$1"\\t"$3;}} }}' | tr -d ">" > {output.bed}
"""

rule MaskFile:
    input:
        bed="wm_mask_intervals.bed",
        asm=assembly
    output:
        masked="assembly.masked.fasta"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/bemask {input.asm} {input.bed} {output.masked}
"""



