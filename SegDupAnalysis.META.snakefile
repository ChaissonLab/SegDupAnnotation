import os
import tempfile
import subprocess
import os.path


# Config
configfile: "/project/mchaisso_100/projects/HPRC/sd_analysis.json"



tempp=config['temp']
if "temp2" not in config:
    config["temp2"] = config["temp"]
    
if config['temp2']!="":
    tempp=config['temp2']


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

assembly="assembly.orig.fasta"
geneModel = config["genemodel"].keys()

#RNADatasets = config["RNAseq"].keys()
spliced=[ "multi", "single"]

bamFiles={f.split("/")[-1]: f for f in config["reads_bam"] }

pos=[]

subs=["all", "high_ident"]

localrules: all, AnnotateResolvedTandemDups, GetUniqueGencodeUnresolvedDupGenes,  IntersectGenesWithFullSDList, FullDupToBed12, FullDupToLinks, MakeWMBed, MaskFile, ConvertHMMCopyNumberToCollapsedDuplications, SortSedef, FilterSedef, CountMaskedSedef, RemoveSedefTooMasked, MakeSedefGraph, MakeSedefGraphTable, FilterByGraphClusters, FullDupToBed12, FiltDupToBed12, GetUniqueGencodeUnresolvedDupGenesCN, GetUniqueGencodeUnresolvedDupGenes, GetGencodeMulticopy, GetGencodeMappedInDup, GetSupportedMulticopy,FindResolvedDupliatedGenes, Bed12ToBed6, CombineGenesWithCollapsedDups, CombineDuplicatedGenes, MinimapGeneModelBed, FilterGencodeBed12, FindGenesInResolvedDups, SelectOneIsoform, SplitSplicedAndSingleExon, AnnotateLowCoverageFlanks, UnionMasked,GetNamedFasta, SelectDups, SortDups, GetDepthOverDups, FilterLowDepthDups, GetFullGeneCountTable, AddCollapsedGenes, GetCombinedTable, SelectDupsOneIsoform, GetFinalMerged, DupsPerContig, GetAllMultiGenes, AnnotateHighIdentity, GetTotalMasked, AnnotateResolvedTandemDups, RemoveBams



#import shutil
#onsuccess:
#    shutil.rmtree(".snakemake")
subs=["all", "high_ident"]

if "name_map" not in config:
    config["name_map"] = "NO_OP"

rule all:
    input:
        fai=assembly+".fai",
        bam=config["bam"],
#        wm_db="wmdb",
#        wm_intv="wm_mask_intervals",
#        masked="assembly.masked.fasta",
        sedef="sedef_out/final.bed",
        sedef_sorted="sedef_out/final.sorted.bed",
        sfilt1="sedef_out/final.sorted.low_copy.bed_stage1",
        sfilt="sedef_out/final.sorted.low_copy.bed",
        sedef_pairs_gml="sedef_out/final.sorted.bed.pairs.gml",
        sedef_pairs_tab="sedef_out/final.sorted.bed.pairs.tab",
        comps="sedef_out/final.sorted.bed.pairs.comps",
        oneIsoform                 = expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt", sub=subs),
        gencodeInDupsNotTandem     = expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.not_tandem", sub=subs),
        realignedOneIsoform        = expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.sam", sub=subs),
        b12                        = expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.sam.bed12", sub=subs),
        realignedOneIsoformFull    = expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.full.sam", sub=subs),
        realignedOneIsoformFullBed = expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.full.sam.bed12", sub=subs),
        geneDups=expand("sedef_out/{sub}/resolved_dups_with_genes.bed", sub=["all", "high_ident"]),
      #  resGeneLinks="circos/genes_in_resolved_dups.links.tsv",
      #  resGeneNames="circos/genes_in_resolved_dups.links.names.tsv",
      #  filtSDResGeneLinks="circos_filtsd/genes_in_resolved_dups.links.tsv",
      #  filtSDResGeneNames="circos_filtsd/genes_in_resolved_dups.links.names.tsv",
        combined_gencode="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined",
        comb_with_unique="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map",
        comb_with_depth="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth",
        gene_count="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.gene_count",
        asm_gene_count="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.asm_gene_count",
        gene_count_2column="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.gene_count_multi_single",        
       # plot="circos/circos.png",
#        plotfilt="circos_filtsd/circos.png",
        splitAndSpliced=expand("sedef_out/{sub}/genes_in_resolved_dups.one_isoform.{sp}.bed", sp=spliced, sub=subs),
        alignedIsoforms=expand("identity.{sp}.bed", sp=spliced),
        sedef_high_ident="sedef_out/high_ident/final.sorted.bed.final.filt",
        sedef_filt="sedef_out/all/final.sorted.bed.final.filt",
        sedef_merged=expand("sedef_out/{sub}/final_filt.merged.bed", sub=subs),
        totalmasked=expand("sedef_out/{sub}/total_masked.txt", sub=subs),
        dupPerContig=expand("sedef_out/{cat}/final_filt.by_contig.bed", cat=subs),
#        repmasked="assembly.repeat_masked.fasta",
#        repmaskedOut="assembly.repeat_masked.fasta.out",
        dups="collapsed_duplications.bed",
        genecol="collapsed_dups_with_genes.bed",
        allDupsWithGenes=expand("sedef_out/{sub}/collapsed_and_resolved_dups_with_genes.bed",sub=subs),
        hmmCopyNumber="hmm/copy_number.tsv",
        rnabam=expand("{data}.mapped.bam",data=geneModel),
        rnabambed=expand("{data}.mapped.bam.bed12",data=geneModel),
        gccn="gencode.mapped.bam.bed12.cn",
        gccn1="gencode.mapped.bam.bed12.cn.one_iso",        
        dupsFromMap=expand("{data}.mapped.bam.bed12.fasta.named.mm2.dups",data=geneModel),
        rnabedFasta=expand("{data}.mapped.bam.bed12.fasta",data=geneModel),
        rnabedFastaNamed=expand("{data}.mapped.bam.bed12.fasta.named",data=geneModel),        
        rnabambed6=expand("{data}.mapped.bam.bed6",data=geneModel),
#        sup=expand("{data}.mapped.bam.bed6.rnaseq-sup",data=geneModel),
#        supdup=expand("{data}.mapped.bam.bed12.dups.sup",data=geneModel),
        colrnabed=expand("{data}.mapped.bam.bed12.dups",data=geneModel),
        resdup=expand("{data}.mapped.resolved_dups.bed",data=geneModel),
        gencodeRes="gencode.mapped.multicopy.bed",
#        gencodeResSup="gencode.mapped.multicopy.bed.supported",
        gencodeSummary="gencode.summary",
        duplicationSummary=expand("sedef_out/{sub}/gencode.combined-duplicated-genes.tsv", sub=subs),
        duplicationSummaryPos=expand("sedef_out/{sub}/gencode.combined-with-position.tsv", sub=subs),
        summaryGencodeResSup=expand("sedef_out/{sub}/gencode.resolved-duplications.tsv", sub=subs),
#        remapped_gencode_cov="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.cov",
        gencodeResIDup="gencode.mapped.multicopy.in_duplication.bed",
#        resdupsup=expand("{data}.mapped.resolved_dups.bed.sup",data=geneModel),
        combineMasked="assembly.union_masked.fasta",
#        RNAseq=expand("RNAseq/{dataset}.bam", dataset=list(config["RNAseq"].keys())),
#        RNAseqCov=expand("RNAseq/{dataset}.bam.cov", dataset=list(config["RNAseq"].keys())),
#        combinedCov="RNAseq/combined.bed",
#        IsoSeq=expand("IsoSeq/{dataset}.bam", dataset=list(config["IsoSeq"].keys())),
        counted="sedef_out/counted.tab",
   #     tandem_dups="sedef_out/tandem_dups.bed",
        low_cov_tandem_dups=expand("sedef_out/{sub}/tandem_dups.low_cov.bed",sub=subs),
        asmMask=expand("{asm}.count_masked", asm=["assembly.orig.fasta", "assembly.masked.fasta", "assembly.repeat_masked.fasta", "assembly.union_masked.fasta"]),
        uniqueDupGenes="gencode.mapped.bam.bed12.dups.unique",
        uniqueDupGenesCN="gencode.mapped.bam.bed12.dups.unique.cn",
       # sdDistPdf=config["species"]+".sd_dist.pdf",
      #  post=dynamic("cn3/post_cn3.{p}.bed"),#,p=pos), #lambda wildcards: getPos("cn3_region.txt")),
        d="done.done",


#
# Simple preprocessing, make sure there is an index on the assembly.
#


rule MakeFaiLinkOrig:
    input:
        asm=config['asm']
    output:
        orig=assembly,
        fai=assembly+".fai"
    params:
        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
ln -s {input.asm} ./{output.orig}
samtools faidx {output.orig}
"""


rule IndexGenome:
    input:
        ref=assembly,
    output:
        gli=assembly+".gli"
    params:
        sd=SD,
        grid_opts=config["grid_medium"],
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
        gli=assembly+".gli"
    output:
        aligned="aligned/{base}.bam"
    params:
        sd=SD,
        ref=assembly,
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

#samtools view -h -F 2304 {input.bam} | samtools fastq - | 
"""

rule MergeBams:
    input:
        aln=expand("aligned/{b}.bam", b=bamFiles.keys())
    output:
        bam=config["bam"]        
    params:
        grid_opts=config["grid_medium"]
    resources:
        load=2
    shell:"""
samtools merge {output.bam} {input.aln} -@2
"""

rule IndexBam:
    input:
        bam=config["bam"]
    output:
        bai=config["bam"] + ".bai"
    resources:
        load=2
    params:
        grid_opts=config["grid_medium"]
    shell:"""
samtools index -@2 {input.bam}
"""

##
## The read depth analysis and duplication assembly all need mapped reads
## 
#rule MakeBam:
#    input:
#        reads=expand("{b}.aligned.bam", b=config["bam"])
#    output:
#        bam=config["bam"]
#    params:
#        sd=SD,
#        ref=assembly,
#        temp=config["temp"],
#        grid_opts=config["grid_large"]
#    shell:"""
#if [ ! -e {params.ref}.gli ]; then
#    {params.sd}/LRA/lra index {params.ref}
#fi
#for f in {input.reads}; do 
#samtools view -h -F 2304 $f | samtools fastq - 
#done | {params.sd}/LRA/lra align {params.ref} /dev/stdin -t 24 -p s | samtools sort -T {params.temp}/asm.$$ -@2 -m4G -o {output.bam}
#samtools index -@2 {output.bam}
#"""
#
#
# Some genomes are poorly masked or are not represented in repeat
# masking libraries. Running windowmasker can identify some repeats
# missed by library based analysis.
#

rule MakeWMDB:
    input:
        asm=assembly
    output:
        wm_db="wmdb"
    params:
        grid_opts=config["grid_medium"]
    resources:
        load=2
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
    resources:
        load=2
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
    resources:
        load=1
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
        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
{params.sd}/bemask {input.asm} {input.bed} {output.masked}
"""

#
# Run repeat masker on the assembly. This will be combined with the
# windowmasker to generate a masked genome.
#

rule SplitGenome:
    input:
        asm="assembly.orig.fasta"
    output:
        split=dynamic("split/to_mask.{region}.fasta")
    params:
        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
mkdir -p split
{params.sd}/DivideFasta.py {input.asm} 10000000 10000 split/to_mask
"""


rule MaskFasta:
    input:
        part="split/to_mask.{region}.fasta"
    output:
        mask="split/to_mask.{region}.fasta.masked"
    params:
        grid_opts=config["grid_medium"],
        repeatLibrary=config["repeat_library"],
        tmpdir=tempp,
    resources:
        load=8
    shell:"""
TEMP="{params.tmpdir}/$$_$RANDOM/"
mkdir -p $TEMP
cp \"{input.part}\" \"$TEMP/to_mask.{wildcards.region}.fasta\" && \
pushd $TEMP &&  \
RepeatMasker {params.repeatLibrary} -pa 8  -s -xsmall \"to_mask.{wildcards.region}.fasta\" && \
popd && \
cp $TEMP/to_mask.\"{wildcards.region}\".fasta.* split/ || true
if [ ! -e split/to_mask.\"{wildcards.region}\".fasta.masked ]; then
  cp split/to_mask.\"{wildcards.region}\".fasta split/to_mask.\"{wildcards.region}\".fasta.masked
fi
rm -rf $TEMP
"""

        
rule MakeToCombine:
    input:
        mask=dynamic("split/to_mask.{region}.fasta.masked"),
    output:
        tc="to_combine.txt"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    run:
        f=open(output.tc,'w')
        f.write("\n".join(input.mask) + "\n")
        f.close()

    
rule CombineMasked:
    input:
        mask="to_combine.txt",
        asm=assembly
    output: 
        masked="assembly.repeat_masked.fasta",
        maskedout="assembly.repeat_masked.fasta.out"
    params:
        grid_opts=config["grid_large"],
        sd=SD
    resources:
        load=1
    shell:"""
{params.sd}/CombineMasked.py {input.mask} {input.asm}.fai {output.masked} {output.maskedout}
"""


#
#  The final masked genome combines wm and repeatmasker
#

if config["repeat_library"] != "pre_masked":
   
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
        resources:
            load=1
        shell:"""
        {params.sd}/comask {output.comb} {input.orig} {input.wm} {input.rm}
        samtools faidx {output.comb}
        """
else:
    rule UnionMasked1:
        input:
            orig="assembly.orig.fasta"
        output:
            comb="assembly.union_masked.fasta"
        shell:"""
ln -sf {input.orig} {output.comb}
samtools faidx {output.comb}
"""
#
# Use excess depth to count duplications
#

rule RunDepthHmm:
    input:
        bam=config["bam"],
        asm="assembly.orig.fasta"
    output:
        vo="hmm/copy_number.tsv",
        cb="hmm/coverage.bins.bed.gz",
        mc="hmm/mean_cov.txt"
    params:
        grid_opts=config["grid_large"],
        sd=SD,
    resources:
        load=16
    shell:"""
snakemake --nolock -p -s {params.sd}/hmm_caller.vert.snakefile -j 16 --rerun-incomplete 
"""

rule RunRefDepthHmm:
    input:
        vo="hmm/copy_number.tsv",
    output:
        vo="hmm_ref/copy_number.tsv",
        cb="hmm_ref/coverage.bins.bed.gz",
        mc="hmm_ref/mean_cov.txt"
    params:
        grid_opts=config["grid_large"],
        sd=SD,
    resources:
        load=16
    shell:"""
snakemake --nolock -p -s {params.sd}/ref_hmm.snakefile -j 16 --rerun-incomplete 
"""

rule ConvertHMMCopyNumberToCollapsedDuplications:
    input:
        bed="hmm/copy_number.tsv"
    output:
        dups="collapsed_duplications.bed"
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
cat {input.bed} | awk '{{ if ($5 > 2) print;}}' > {output.dups}
"""

rule MakeCoverageBins:
    input:
        cb="collapsed_duplications.bed"
    output:
        cbcol="collapsed_duplications.bed.collapse",
        s="pre.collapsed_duplications.split.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
bedtools merge -i {input.cb} -c 5 -o collapse > {output.cbcol}
{params.sd}/SplitCoverageBins.py {output.cbcol} | bedtools merge -c 4,4 -o min,max > {output.s}
"""

rule Postcn3:
    input:
        s="pre.collapsed_duplications.split.bed"
    output:
        pre="pre_cn3.txt",
        reg="cn3_region.txt",
        nf="cn3.nucfreq.bed.gz"
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
        nf="cn3.nucfreq.bed.gz",
        reg="cn3_region.txt",
     #   done="getPos.done",
       # ps=lambda wildcards: pos[wildcards.p],
    output:
        post="cn3/post_cn3.bed",
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
        post="cn3/post_cn3.bed",
        #expand("cn3/post_cn3.{p}.bed",p=lambda wildcards: getPos("cn3_region.txt")),
        s="pre.collapsed_duplications.split.bed"
    output:
        ss="collapsed_duplications.split.bed"
    resources:
        load=1
    params:
        grid_opts=config["grid_small"],
    shell:"""
intersectBed -v -a {input.s} -b <( cat {input.post} |grep fail) > {output.ss} 
    """





#
# The following rules do the initial resolved repeat detection with
# sedef, and then postprocess the output to remove excess duplications.
# 

#
# Initial run of sedef.
#
rule RunSedef:
    input:
        asm="assembly.union_masked.fasta"
    output:
        done="sedef_out/final.bed"
    params:
        grid_opts=config["grid_sedef"],
        sd=SD
    resources:
        load=12
    shell:"""

module load gcc/8.3.0
module load time
module load parallel
export PATH=$PATH:{params.sd}/sedef

{params.sd}/sedef.sh  {input.asm} -j 12
"""
#
# Sort by chrom and start, fixing a bug in sedef output that misses
# column output.
#

rule SortSedef:
    input:
        bed="sedef_out/final.bed"
    output:
        s="sedef_out/final.sorted.bed"
    params:
        grid_opts=config["grid_medium"]
    resources:
        load=1
    shell:"""
first=`head -1 {input.bed} | awk '{{ print NF;}}'`
sort -k1,1 -k2,2n {input.bed} | \
  awk -vf=$first '{{ if (NF == f) print;}}' | \
  bedtools groupby -g 1-6 -o first -full -c 1 | cut -f 1-$first > {output.s}
"""

#
# Sedef annotates identity of repeat, filter that here.
#

rule FilterSedef:
    input:
        bed="sedef_out/final.sorted.bed"
    output:
        filtered="sedef_out/final.sorted.first_filtered.bed",
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
awk  '{{ if ($NF > 0.5) {{ print;}} }} ' > {output.filtered} < {input.bed}
"""

#
# Some of the sedef alignments are masked repeats. Remove those.
#
rule CountMaskedSedef:
    input:
        s="sedef_out/final.sorted.first_filtered.bed",
        asm="assembly.union_masked.fasta",
    output:
        f="sedef_out/final.sorted.bed.frac_masked"
    params:
        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
{params.sd}/CountMaskedBed.py {input.s} {input.asm} > {output.f}
"""

#
# Count how many times an entry appears 
#

#
# Count of often an entry overlaps with others
#

rule CountRepeatedEntries:
    input:
        missed="sedef_out/final.sorted.first_filtered.bed"
    output:
        counted="sedef_out/counted.tab"
    params:
        grid_opts=config["grid_medium"]
    resources:
        load=1
    shell:"""
bedtools intersect -loj -a {input.missed} -b {input.missed} -sorted -f 0.9 -r | bedtools groupby -g 1-6 -c 2 -o count | awk '{{ print $NF;}}' > {output.counted}
"""


#
# Remove too masked, or too high copy
#
rule RemoveSedefTooMasked:
    input:
        counted="sedef_out/counted.tab",
        fracMasked="sedef_out/final.sorted.bed.frac_masked",
        inputBed="sedef_out/final.sorted.first_filtered.bed"
    output:
        filt1="sedef_out/final.sorted.low_copy.bed_stage1",
        filt="sedef_out/final.sorted.low_copy.bed"
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
paste  {input.inputBed} {input.counted} {input.fracMasked} | awk '{{ if  ($NF < 0.9 ) print;}}'  > {output.filt1}
paste  {input.inputBed} {input.counted} {input.fracMasked} | awk '{{ if ($(NF-1) < 20 && $NF < 0.9 ) print;}}'  > {output.filt}
"""

#
# Use graph analysis to finally remove high copy dups.
#
rule MakeSedefGraph:
    input:
        pairs="sedef_out/final.sorted.low_copy.bed"
    output:
        gml="sedef_out/final.sorted.bed.pairs.gml",
        tab="sedef_out/final.sorted.bed.pairs.tab",  
    params:
        grid_opts=config["grid_medium"],
        sd=SD
    resources:
        load=1
    shell:"""
cat {input.pairs} | grep -v "^#" | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\tLEFT"; print $4"\\t"$5"\\t"$6"\\t"$1"\\t"$2"\\t"$3"\\tRIGHT";}}' | bedtools sort > {input.pairs}.lr
{params.sd}/CountResolvedDuplicationMultiplicity.py {input.pairs}.lr {output.gml} {output.tab} --overlap 0.7


"""

#
# Transform left/right (source/dest clusters) sedef output to the same # rows as input
#
rule MakeSedefGraphTable:
    input:
        pairs="sedef_out/final.sorted.low_copy.bed",
        tab="sedef_out/final.sorted.bed.pairs.tab"     
    output:
        comps="sedef_out/final.sorted.bed.pairs.comps"
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
paste {input.pairs}.lr {input.tab} | grep LEFT >  {output.comps}
"""
#
#
# Filter based on graph
#
rule FilterByGraphClusters:
    input:
        bed="sedef_out/final.sorted.low_copy.bed",
        comps="sedef_out/final.sorted.bed.pairs.comps"
    output:
        filt="sedef_out/all/final.sorted.bed.final.filt",
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
mkdir -p sedef_out/all
paste {input.bed} {input.comps} | bioawk -c hdr '{{ if ($(NF-1) <= 20) print;}}' > {output.filt}

"""

rule AnnotateHighIdentity:
    input:
        filt="sedef_out/all/final.sorted.bed.final.filt",
    output:
        sedef_high_ident="sedef_out/high_ident/final.sorted.bed.final.filt",
    shell:"""
mkdir -p sedef_out/high_ident
cat {input.filt} | awk '{{ if ($21 >= 0.98) print;}}' > {output.sedef_high_ident}
"""

rule GetTotalMasked:
    input:
        filt="sedef_out/{cat}/final.sorted.bed.final.filt",
    output:
        total="sedef_out/{cat}/total_masked.txt"
    shell:"""
cat {input.filt} | awk '{{ print $1"\\t"$2"\\t"$3"\\t"; print $4"\\t"$5"\\t"$6; }}'  | bedtools sort | bedtools merge | awk '{{ s+= $3-$2;}} END {{ print s;}}' > {output.total}
"""

rule GetFinalMerged:
    input:
        filt="sedef_out/{cat}/final.sorted.bed.final.filt",
    output:
        merged="sedef_out/{cat}/final_filt.merged.bed"
    shell:"""
cat {input.filt} | awk '{{ print $1"\\t"$2"\\t"$3; print $4"\\t"$5"\\t"$6;}}' | bedtools sort | bedtools merge > {output.merged}
"""

rule DupsPerContig:
    input:
        merged="sedef_out/{cat}/final_filt.merged.bed",
        asm="assembly.orig.fasta"        
    output:
        dupPerContig="sedef_out/{cat}/final_filt.by_contig.bed",
    shell:"""
set +e
while read line; do
    name=`echo $line | awk '{{ print $1;}}' `
    len=`echo $line | awk '{{ print $2; }}' `
    grep -q $name {input.merged}
    if [ $? -eq 0 ]; then
      n=`grep $name {input.merged} | awk '{{ s+=$3-$2; }} END {{ print s;}}'`
    else
      n="0"
    fi
    echo -e "$name\\t$len\\t$n"
done < {input.asm}.fai > {output.dupPerContig}
"""


rule AnnotateResolvedTandemDups:
    input:
        repeats="sedef_out/{seg}/final.sorted.bed.final.filt",
    output:
        tandem="sedef_out/{seg}/tandem_dups.bed"

    shell:"""
cat  {input.repeats} | awk '{{ if ($1 == $4 && $5 > $3 && $5-$3 < 1000 && prevPos!=$2) {{ print;}}  prevPos=$2;}}'  > {output.tandem}

"""


rule AnnotateLowCoverageFlanks:
    input:
        tandem="sedef_out/{seg}/tandem_dups.bed",
        bam=config["bam"],
        bai=config["bam"] + ".bai",
        meancov="hmm/mean_cov.txt"
    output:
        low_cov_tandem_dups="sedef_out/{seg}/tandem_dups.low_cov.bed",
    params:
        grid_opts=config["grid_medium"]
    resources:
        load=1
    shell:"""
cut -f 1-3 {input.tandem} > {input.tandem}.3
minCov=`cat {input.meancov} | awk '{{ print $1/5}}'`
while read -r line; do 
    region=`echo $line | awk '{{ print $1":"$3-100"-"$3+100; }}'`
    minCov=`samtools depth {input.bam} -r $region | awk 'BEGIN{{ m=99999999;}} {{ if ($3 < m) {{ m=$3;}} }} END {{ print m;}}'`
    echo -e "$line\t$minCov" | tr " " "\t"
done < {input.tandem}.3  | awk -v m=$minCov '{{ if ($4 < m) print;}}' > {output.low_cov_tandem_dups}
"""


#END
#rule MakeSymmetricalDups:
#    input:
#        dups="sedef_out/{sub}/final.sorted.bed.final.filt",
#    output:
#        dups="sedef_out/{sub}/final.sorted.bed.final.filt.sym",
#    shell:"""
#cat {input.dups} | awk '{{ 
#""

rule GencodeCN:
    input:
        gc="gencode.mapped.bam.bed12",
        cn="hmm/copy_number.tsv",
        genome="assembly.orig.fasta"
    params:
        grid_opts=config["grid_small"],
    output:
        gccn="gencode.mapped.bam.bed12.cn",
    shell:"""
bedtools intersect -loj -g {input.genome}.fai -a gencode.mapped.bam.bed12 -b hmm/copy_number.tsv -sorted | awk '{{ if (NF== 17) print;}}' | bedtools groupby -g 1-4 -c 17 -o mean > {output.gccn}

"""

rule GencodeOneIsoCN:
    input:
        gccn="gencode.mapped.bam.bed12.cn",
    output:
        gccn_iso="gencode.mapped.bam.bed12.cn.one_iso",
    params:
        grid_opts=config["grid_small"],
    shell:"""
cat {input.gccn} | tr "|" "\\t" | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$9"\\t"$NF;}}' | bedtools groupby -g 4 -c 5 -o median > {output.gccn_iso}
"""

rule IntersectGenesWithFullSDList:
    input:
        dups="sedef_out/{sub}/final.sorted.bed.final.filt",
        refGenes="gencode.mapped.bam.bed12"
    output:
        fullDup="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.full.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
na=`head -1 {input.dups} | awk '{{ print NF;}}'`
nb=`head -1 {input.refGenes} | awk '{{ print NF;}}'`
tot=`echo "" | awk -va=$na -vb=$nb '{{ print a+b;}}'`
bedtools intersect -a {input.refGenes} -b {input.dups} -f 1 -wb | \
 awk -vt=$tot '{{ if (NF == t) print;}}' | \
 bedtools groupby -g 1-4 -c 1 -o first -full | cut -f 1-$tot | \
 {params.sd}/FilterMembersFromSameIsoformSet.py stdin | \
 {params.sd}/FilterLongestInOverlapSet.py stdin > {output.fullDup}
"""



rule ResolvedDupsWithGenes:
    input:
        dupGenes="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.full.bed",
        dups="sedef_out/{sub}/final.sorted.bed.final.filt",
    output:
        geneDups="sedef_out/{sub}/resolved_dups_with_genes.bed"
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -a {input.dups} -b {input.dupGenes} -F 1 -wb | \
   awk 'BEGIN {{OFS="\\t";}} {{ print $1,$2,$3,$42,"resolved"; print $4,$5,$6,$42,"resolved"}}' | \
   bedtools sort |  bedtools merge -c 4,5 -o collapse >  {output.geneDups}
"""

rule CombineCoordinatesOfDupsWithGenes:
    input:
        collGene="collapsed_dups_with_genes.bed",
        resGene="sedef_out/{sub}/resolved_dups_with_genes.bed"
    output:
        allDupsWithGenes="sedef_out/{sub}/collapsed_and_resolved_dups_with_genes.bed",
    params:
        grid_opts=config["grid_small"]
    shell:"""
cat {input.resGene} {input.collGene} > {output.allDupsWithGenes}
"""


rule SummarizeGencode:
    input:
        bam="gencode.mapped.bam",
        bed="gencode.mapped.bam.bed6",
        bed12="gencode.mapped.bam.bed12",
    output:
        res="gencode.summary"
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:"""
nl=`samtools view {input.bam} | wc -l | awk '{{ print $1;}}'`;
nb=`cat {input.bed} | awk '{{ n+=$3-$2;}} END {{ print n;}}'`;
ng=`cut -f 4 {input.bed12} | sort | uniq | awk 'END {{ print NR;}}'`
echo -e "Lines\t$nl" > {output.res}
echo -e "Bases\t$nb" >> {output.res}
echo -e "Genes\t$ng" >> {output.res}
"""

rule RealignFullDup:
    input:
        oneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.full.bed",
        asm="assembly.orig.fasta",
        refGenes=config["genemodel"]["gencode"],
        refGeneBam="gencode.mapped.bam"
    output:
        realignedOneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.full.sam",
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
{params.sd}/RemapRegion.sh {input.oneIsoform} {input.asm} {input.refGenes} {input.refGeneBam} {output.realignedOneIsoform}
"""

rule FullDupToBed12:
    input:
        realignedOneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.full.sam",
    output:
        bed="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.full.sam.bed12",
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
bedtools bamtobed -bed12 -i {input.realignedOneIsoform} | {params.sd}/SimplifyNameInBed12.py | sort -k4,4 -k1,1 -k2,2n | bedtools groupby -g 1-3 -c 2 -o first -full > {output.bed}
"""

rule FiltDupToBed12:
    input:
        realignedOneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.full.sam",
    output:
        bed="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.sam.bed12",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
bedtools bamtobed -bed12 -i {input.realignedOneIsoform} | {params.sd}/SimplifyNameInBed12.py | sort -k4,4 -k1,1 -k2,2n | bedtools groupby -g 1-3 -c 2 -o first -full > {output.bed}
"""

#rule FiltDupToLinks:
#    input:
#        bed="sedef_out/high_ident/genes_in_resolved_dups.one_isoform.bed.filt.sam.high_ident.bed12",
#    output:
#        links="circos_filtsd/genes_in_resolved_dups.links.tsv",
#        names="circos_filtsd/genes_in_resolved_dups.links.names.tsv",
#    params:
#        sd=SD,
#        grid_opts=config["grid_small"]
#    resources:
#        load=1
#    shell:"""
#mkdir -p circos
#cat {input.bed} | tr "#" "_" | {params.sd}/Bed12ToArcs.py {output.links} {output.names}
#"""


#rule FullDupToLinks:
#    input:
#        bed="sedef_out/high_ident/genes_in_resolved_dups.one_isoform.bed.full.sam.bed12",
#    output:
#        links="circos/genes_in_resolved_dups.links.tsv",
#        names="circos/genes_in_resolved_dups.links.names.tsv"
#    params:
#        sd=SD,
#        grid_opts=config["grid_small"]
#    resources:
#        load=1
#    shell:"""
#mkdir -p circos
#cat {input.bed} | tr "#" "_" |  {params.sd}/Bed12ToArcs.py {output.links} {output.names}
#"""

rule RemapBed:
    input:
        oneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt",
        asm="assembly.orig.fasta",
        refGenes=config["genemodel"]["gencode"],
        refGeneBam="gencode.mapped.bam"
    output:
        realignedOneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.sam",
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
{params.sd}/RemapRegion.sh {input.oneIsoform} {input.asm} {input.refGenes} {input.refGeneBam} {output.realignedOneIsoform}
"""
   
rule RemappedBedToSummaryTable:
    input:
        realignedOneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.sam",
    output:
        b12="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.sam.filt.bed12",
        counts="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.sam.filt.counts",
    params:
        sd=SD,
        grid_opts=config["grid_medium"]
    resources:
        load=1
    shell:"""
samtools view {input.realignedOneIsoform} -o {input.realignedOneIsoform}.bam
bedtools bamtobed -bed12 -i {input.realignedOneIsoform}.bam | \
  {params.sd}/SimplifyNameInBed12.py | \
  sort -k4,4 -k1,1 -k2,2n | \
  bedtools groupby -g 1-3 -c 2 -o first -full | {params.sd}/RemoveOverlappingGenes.py >  {output.b12}


grep -v "^@" {input.realignedOneIsoform} | awk '{{ for (i=1; i <= NF; i++) {{ if (substr($i,0,3) == "GR:") {{ print substr($i,6);}} }} }}' > groups.txt

bedtools groupby -g 4 -c 4 -i {output.b12} -o count > {output.counts}
"""

#rule MakeCircOS:
#    input:
#        asm="assembly.orig.fasta",
#       coll="gencode.mapped.bam.bed12.dups.unique",
#        links="circos/genes_in_resolved_dups.links.tsv",
#        names="circos/genes_in_resolved_dups.links.names.tsv"
#    output:
#        plt="circos/circos.png"
#    params:
#        sd=SD,
#        grid_opts=config["grid_medium"],
#        name_map=config["name_map"]
#    resources:
#       load=1
#    shell:"""
#mkdir -p circos
#cat {input.coll} | tr "#" "_" | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$(NF-1)"\\t"$NF;}}' > {input.coll}.name.cn
#cat {input.asm}.fai | tr "#" "_" > {input.asm}.fai.u
#cut -f 1 {input.asm}.fai.u | sed "s/_/__/g" | tr "#" "_" > {input.asm}.display_name
#paste {input.asm}.fai.u {input.asm}.display_name | awk '{{ if ($2 > 200000) {{ print "chr\\t-\\tvar"$1"\\t"$6"\\t"0"\\t"$2"\\t"$1;}}  }}' | \
 # {params.sd}/ReplaceChromName.py {params.name_map} >  circos/karyotype.txt

#cat {input.coll}.name.cn | awk '{{ print "var"$1"\\t"$2"\\t"$3"\\t"$4;}}' | {params.sd}/ReplaceChromName.py {params.name_map} > circos/cn.lab.txt
#cat {input.coll}.name.cn | awk '{{ print "var"$1"\\t"$2"\\t"$3"\\t"$5;}}' | {params.sd}/ReplaceChromName.py {params.name_map} > circos/cn.txt

#{params.sd}/MakeDup.py --bed {input.coll}  --collapsed {input.coll}.name.cn --links circos/resolved_dups.txt --labels circos/resolved_dups.labels.txt
#cd circos && circos --conf {params.sd}/circos.conf 
#"""


#rule MakeCircOSHighIdentityDups:
#    input:
#        asm="assembly.orig.fasta",
#        coll="gencode.mapped.bam.bed12.dups.unique",
#        links="circos_filtsd/genes_in_resolved_dups.links.tsv",
#        names="circos_filtsd/genes_in_resolved_dups.links.names.tsv"
#    output:
#        plt="circos_filtsd/circos.png"
#    params:
#        sd=SD,
#        grid_opts=config["grid_medium"]
#    resources:
#        load=1
#    shell:"""
#mkdir -p circos_filtsd
#
#cat {input.coll} | tr "#" "_" | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$(NF-1)"\\t"$NF;}}' > {input.coll}.name.cn
#cat {input.asm}.fai | tr "#" "_" > {input.asm}.fai.u
#cut -f 1 {input.asm}.fai.u | sed "s/_/__/g" > {input.asm}.display_name
#paste {input.asm}.fai.u {input.asm}.display_name | awk '{{ if ($2 > 200000) {{ print "chr\\t-\\tvar"$1"\\t"$6"\\t"0"\\t"$2"\\t"$1;}} }}' > circos_filtsd/karyotype.txt
#
#cat {input.coll}.name.cn | awk '{{ print "var"$1"\\t"$2"\\t"$3"\\t"$4;}}' > circos_filtsd/cn.lab.txt
#cat {input.coll}.name.cn | awk '{{ print "var"$1"\\t"$2"\\t"$3"\\t"$5;}}' > circos_filtsd/cn.txt
#
#{params.sd}/MakeDup.py --bed  --collapsed {input.coll}.name.cn --links circos_filtsd/resolved_dups.txt --labels circos_filtsd/resolved_dups.labels.txt
#cd circos_filtsd && circos --conf {params.sd}/circos.conf 
#
#
#"""

rule GetAllMultiGenes:
    input:
        single="gencode.mapped.bam.bed12.dups.unique.cn",
        multi="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined"
    output:
        both="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map"
    params:
        sd=SD
    shell:"""
{params.sd}/AppendCollapsedDupsList.py  {input.multi} {input.single} | sort -k4,4 > {output.both}
"""

rule GetUniqueGencodeUnresolvedDupGenesCN:
    input:
        uniqueDupGenes="gencode.mapped.bam.bed12.dups.unique"
    output:
        uniqueDupGenesCN="gencode.mapped.bam.bed12.dups.unique.cn"
    params:
        grid_opts=config["grid_small"],
    resources:
        load=1
    shell:"""
cut -f 1-4,16 {input.uniqueDupGenes} > {output.uniqueDupGenesCN}
"""
    
rule GetUniqueGencodeUnresolvedDupGenes:
    input:
        bed="gencode.mapped.bam.bed12.dups"
    output:
        unique="gencode.mapped.bam.bed12.dups.unique"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
cat {input.bed} |  {params.sd}/WriteNonOverlapping.py | {params.sd}/SimplifyNameInBed12.py | sort -k4,4 -k1,1 -k2,2n > {output.unique}
"""
   

rule GetGencodeMulticopy:
    input:
        bed="gencode.mapped.bam.bed12.iso_filt",
    output:
        gencodeRes="gencode.mapped.multicopy.bed"
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    resources:
        load=1
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
        dups="sedef_out/high_ident/final.sorted.bed.final.filt"
    output:
        inDup="gencode.mapped.multicopy.in_duplication.bed"
    params:
        grid_opts=config["grid_small"],
        load=1
    shell:"""
cat {input.dups} | awk '{{ if ($8 <= 10) print;}}' > {input.dups}.ten
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
    resources:
        load=1
    shell:"""
bedtools intersect -u -a {input.bed}  -b {input.sup} > {output.sup}
"""


rule SummarizeMulticopy:
    input:
        sup="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.sam.filt.bed12",
    output:
        summary="sedef_out/{sub}/gencode.resolved-duplications.tsv",
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
cat {input.sup} | \
  awk '{{ print $13"\\t"$14"\\t"$15"\\t"$4; print $16"\\t"$17"\\t"$18"\\t"$4;}}' | \
  sort -k4,4 -k1,1 -k2,2n  | \
  bedtools groupby -g 1-4 -o first -full -c 4 | cut -f 4 | uniq -c | awk '{{ print $2"\\t"$1"\\tresolved"; }}' > {output.summary}
"""



#
# Find all the full-length genes that overlap a resolved duplication.
#
rule FindResolvedDupliatedGenes:
    input:
        rnabed="{data}.mapped.bam.bed12.iso_filt",
        sedef="sedef_out/high_ident/final.sorted.bed.final.filt",
    output:
        resdup="{data}.mapped.resolved_dups.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""

na=`head -1 {input.rnabed} | awk '{{ print NF;}}'`
nb=`head -1 {input.sedef} | awk '{{ print NF;}}'`
tot=`echo "" | awk -va=$na -vb=$nb '{{ print a+b;}}'`
bedtools intersect -a {input.rnabed} -b {input.sedef} -loj -f 1 | \
  awk -vt=$tot '{{ if (NF == t) print; }}' | \
  awk '{{ if ($13 != ".") print; }}' | \
  bedtools sort >  {output.resdup}
"""



rule Bed12ToBed6:
    input:
        pre="{base}.bed12"
    output:
        post="{base}.bed6"
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
bedtools bed12tobed6 -i {input.pre} > {output.post}
"""

#rule GetRNASeqSupportedTranscript:
#    input:
#        bed6="{base}.bed6",
#        comb="RNAseq/combined.bed"
#    output:
#        sup="{base}.bed6.rnaseq-sup"
#    params:
#        grid_opts=config["grid_small"]
#    shell:"""
#bedtools intersect -a {input.bed6} -b {input.comb} -wao | bedtools sort > {output.sup}
#"""
#
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
    resources:
        load=1
    shell:"""
if [ ! -e {input.asm}.fai ]; then
 samtools faidx {input.asm}
fi
na=`head -1 {input.dups} | awk '{{ print NF;}}'`
nb=`head -1 {input.rnabed} | awk '{{ print NF;}}'`
tot=`echo "" | awk -va=$na -vb=$nb '{{ print a+b;}}'`

bedtools intersect -f 1 -g {input.asm}.fai -sorted -loj -a {input.rnabed} -b {input.dups} | \
awk -vt=$tot '{{ if (NF == t) print; }}' | \
 awk '{{ if ($13 != ".") print;}}'| \
 bedtools sort | \
 {params.sd}/FilterMembersFromSameIsoformSet.py stdin | \
  {params.sd}/FilterLongestInOverlapSet.py stdin > {output.rnabedout}

"""

#
# Once genes have been identified in collapses, this annotates which collapses have genes.
#
rule GetCollapseWithGenes:
    input:
        genes="gencode.mapped.bam.bed12.dups",
        dups="collapsed_duplications.split.bed",
    output:
        genecol="collapsed_dups_with_genes.bed",
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:"""
bedtools intersect -a {input.dups} -b {input.genes} -F 1 -wb | \
   awk '{{ print $1"\\t"$2"\\t"$3"\\t"$9"\\tcollapsed/"$4;}}' | \
  {params.sd}/SimplifyName.py > {output.genecol}
"""


rule SummarizeGenesWithCollapsedDups:
    input:
        genes="gencode.mapped.bam.bed12.dups.unique",
    output:
        coll="gencode.collapsed.tsv"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
cat {input.genes} | \
cut -f 4,16 |  awk '{{ print $1"\\t"$2"\\tcollapse";}}' > {output.coll}
"""

rule CombineDuplicatedGenesWithCoordinates:
    input:
        col = "gencode.mapped.bam.bed12.dups",
        res = "sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.sam.filt.bed12"
    output:
        comb= "sedef_out/{sub}/gencode.combined-with-position.tsv"
    params:
        sd=SD, 
        grid_opts=config["grid_small"]
    shell:"""
cat gencode.mapped.bam.bed12.dups |  \
  {params.sd}/SimplifyNameInBed12.py | \
  awk '{{ print $4"\\t"$(NF-1) "\\t" $1 "\\t" $2 "\\t" $3"\\t1\\tcollapse";}}' > {output.comb}.tmp

cat {input.res} | \
   awk 'BEGIN {{p="";}} {{ if ($4 != p) {{ i=1; }} else {{ i+=1;}} print $4 "\\t1\\t" $1 "\\t" $2 "\\t" $3 "\\t" i"\\tresolved"; p=$4;}}' >> {output.comb}.tmp
sort {output.comb}.tmp > {output.comb}
rm {output.comb}.tmp

"""

rule CombineDuplicatedGenes:
    input:
        coll = "gencode.collapsed.tsv",
        res  = "sedef_out/{sub}/gencode.resolved-duplications.tsv",
    output:
        summary = "sedef_out/{sub}/gencode.combined-duplicated-genes.tsv"
    params:
        grid_opts = config["grid_small"]
    resources:
        load=1
    shell:"""
cat {input.coll} {input.res} | sort > {output.summary}
"""

#rule GetSupportedDups:
#    input:
#        dup="{data}.mapped.bam.bed12.dups",
#        rna="{data}.mapped.bam.bed6.rnaseq-sup"
#    output:
#        supdup="{data}.mapped.bam.bed12.dups.sup",
#    params:
#        grid_opts=config["grid_small"]
#    shell:"""
#bedtools intersect -sorted -wao -a {input.dup} -b {input.rna} | bedtools groupby -c 1 -o first -full > {output.supdup}
#"""
#
#rule GetSupportedResolvedDups:
#    input:
#        dup="{data}.mapped.resolved_dups.bed",
#        rna="{data}.mapped.bam.bed6.rnaseq-sup"
#    output:
#        supdup="{data}.mapped.resolved_dups.bed.sup",
#    params:
#        grid_opts=config["grid_small"]
#    shell:"""
#bedtools intersect -sorted -wao -a {input.dup} -b {input.rna} | bedtools groupby -c 1 -o first -full > {output.supdup}
#"""


#rule CountDuplicatedGenes:
#    input:
#        rnabed=
#


rule CountMaskedAsmp:
    input:
        asm="{assembly}"
    output:
        countMasked="{assembly}.count_masked"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
cat {input.asm} | {params.sd}/nl > {output.countMasked}
"""



rule FilterGencodeBed12:
    input:
        gencode="{data}.mapped.bam.bed12"
    output:
        iso_filt="{data}.mapped.bam.bed12.iso_filt"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
if [ {wildcards.data} = "gencode" ]; then
  {params.sd}/FilterMembersFromSameIsoformSet.py {input.gencode} |  {params.sd}/FilterLongestInOverlapSet.py  /dev/stdin   | \
   egrep "CDS|protein_coding" > {output.iso_filt}
else
  {params.sd}/FilterLongestInOverlapSet.py   {input.gencode}  > {output.iso_filt}
fi
"""

rule FindGenesInResolvedDups:
    input:
        notmasked="sedef_out/{sub}/final.sorted.bed.final.filt",
        gencode="gencode.mapped.bam.bed12"
    output:
        gencodeInDups="sedef_out/{sub}/genes_in_resolved_dups.bed.filt",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""

na=`head -1 {input.notmasked} | awk '{{ print NF;}}'`
nb=`head -1 {input.gencode} | awk '{{ print NF;}}'`
tot=`echo "" | awk -va=$na -vb=$nb '{{ print a+b+1;}}'`

cat '{input.notmasked}' | \
  awk '{{ print $0"\\t"NR; a=$1;b=$2;c=$3;d=$4;e=$5;f=$6; $1=d;$2=e;$3=f; $4=a;$5=b;$6=c; print $0"\\t"NR;}}' | \
    tr " " "\\t" | \
    bedtools intersect -loj -a {input.gencode} -b stdin -f 1  | awk '{{ if ($13 != ".") print ;}}'  | awk -vt=$tot '{{ if (NF==t) print;}}' | \
 {params.sd}/FilterMembersFromSameIsoformSet.py stdin | \
  {params.sd}/FilterLongestInOverlapSet.py stdin >  {output.gencodeInDups}

"""

rule SelectOneIsoform:
    input:
        anyIsoform="sedef_out/{sub}/genes_in_resolved_dups.bed.filt",
        falseTandem="sedef_out/{sub}/tandem_dups.low_cov.bed",
    output:
        oneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt",
        gencodeInDupsNotTandem     = "sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt.not_tandem",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
na=`head -1 {input.anyIsoform} | awk '{{ print NF;}}'`
nb=`head -1 {input.falseTandem} | awk '{{ print NF;}}'`
tot=`echo "" | awk -va=$na -vb=$nb '{{ print a+b;}}'`

{params.sd}/FilterMembersFromSameIsoformSet.py {input.anyIsoform}  > {output.oneIsoform}
bedtools intersect -v -a {output.oneIsoform} -b  {input.falseTandem} > {output.gencodeInDupsNotTandem}
"""

rule SplitSplicedAndSingleExon:
    input:
        oneIsoform="sedef_out/{sub}/genes_in_resolved_dups.one_isoform.bed.filt",
    output:
        multiAndOne=expand("sedef_out/{{sub}}/genes_in_resolved_dups.one_isoform.{sp}.bed", sp=spliced),
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
{params.sd}/SplitSplicedGenes.py {input.oneIsoform} {output.multiAndOne}
"""
    
rule AlignIsoforms:
    input:
        asm="assembly.union_masked.fasta",
        refGenes=config["genemodel"]["gencode"],
        dups="sedef_out/high_ident/genes_in_resolved_dups.one_isoform.{sp}.bed"
    output:
        alignedIsoforms="identity.{sp}.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
{params.sd}/AlignIsoforms.sh {input.dups} {input.refGenes} {input.asm} > {output.alignedIsoforms}
"""



rule ExtractRepeatsFromCounted:
    input:
        counted="sedef_out/counted.bed"
    output:
        reps="sedef_out/counted.bed.rep"
    params:
        grid_opts=config["grid_small"]
    resources:
        load=1
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
    resources:
        load=1
    shell:"""
cat {input.reps} | awk '{{ print $1":"$2"-"$3;}}' > {output.rgn}
"""

rule CountRepeatIdentity:
    input:
        reps="sedef_out/counted.bed.rep.rgn",
        asm="assembly.union_masked.fasta"
    output:
        ident="sedef_out/counted.bed.rep.rgn.ident"
    params:
        grid_opts=config["grid_small"],
        sd=SD,    
    resources:
        load=1
    shell:"""
samtools faidx {input.asm} -r {input.reps} | {params.sd}/nl > {output.ident}
"""
    
#rule MapRNASeq:
#    input:
#        rna=lambda wildcards : config["RNAseq"][wildcards.dataset],
#        asm="assembly.orig.fasta"
#    output:
#        aln="RNAseq/{dataset}.bam"
#    params:
#        grid_opts=config["grid_large"]
#    shell:"""
#mkdir -p RNAseq
#minimap2 -t 16 -ax splice -uf -C5 {input.asm} {input.rna} | samtools view -uS - | samtools sort -T $TMPDIR/{wildcards.dataset}.$$ -m2G -o {output.aln}
#samtools index {output.aln}
#"""

#rule RNASeqCov:
#    input:
#        aln="RNAseq/{dataset}.bam",
#        asm="assembly.orig.fasta"
#    output: 
#        cov="RNAseq/{dataset}.bam.cov"
#    params:
#        grid_opts=config["grid_large"],
#        sd=SD
#    shell:"""
#cat {input.asm}.fai | awk '{{ if ($2 > 20000) {{ print $1":1-"$2;}} }}'> {input.asm}.rgn
#
#cat {input.asm}.rgn | xargs -P16 -I {{}} {params.sd}/BamToFreq.sh {input.aln} {{}} {input.asm} {output.cov}
#cat {output.cov}.*.rgn | bedtools sort > {output.cov} 
#
#"""    

#rule CombineRNASeqCov:
#    input:
#        rna=expand("RNAseq/{dataset}.bam.cov", dataset=config["RNAseq"].keys())
#    output:
#        combinedCov="RNAseq/combined.bed",
#    params:
#        grid_opts=config["grid_medium"]
#    shell:"""
#cat {input.rna} | bedtools sort | bedtools merge -c 4 -o max > {output.combinedCov}
#"""
#
rule MapIsoSeq:
    input:
        rna=lambda wildcards : config["IsoSeq"][wildcards.dataset],
        assembly="assembly.orig.fasta"
    output:
        aln="IsoSeq/{dataset}.bam" 
    params:
        grid_opts=config["grid_large"],
        temp=config["temp"],
    resources:
        load=16
    shell:"""
mkdir -p IsoSeq
minimap2 -t 16 -ax splice -uf -C5 {input.assembly} {input.rna} | samtools view -uS - | samtools sort -T {params.temp} -m2G -o {output.aln}
samtools index -c {output.aln}
"""

    
rule FilterPSL:
    input:
        psl="{base}.psl"
    output:
        filt="{base}.psl.filt"
    params:
        grid_opts=config["grid_blat"],
        sd=SD
    resources:
        load=1
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
    resources:
        load=1
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
        sd=SD,
        temp=config["temp"],
    resources:
        load=8
    shell:"""
minimap2 -x splice -a -t 4 {input.asm} {input.fa}  | samtools sort -T {params.temp}/tmp.$$ -o {output.bam}
samtools index -c {output.bam}
"""

rule FilterGeneModel:
    input:
        mmbam="{data}.minimapped.bam"
    output:
        bam="{data}.mapped.bam"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    resources:
        load=1
    shell:"""
{params.sd}/FilterMappedLength.py {input.mmbam} any | samtools view -b -o {output.bam}
samtools index -c {output.bam}
"""

rule GetGeneBoundaryFasta:
    input:
        bed="{data}.mapped.bam.bed12",
        asm="assembly.orig.fasta"
    output:
        fa="{data}.mapped.bam.bed12.fasta",        
    params:
        grid_opts=config["grid_small"],
    shell:"""
cat {input.bed} | awk '{{ print $1":"$2"-"$3;}}' > {input.bed}.rgn
samtools faidx {input.asm} -r {input.bed}.rgn > {output.fa}
"""

rule GetNamedFasta:
    input:
        fa="{data}.mapped.bam.bed12.fasta",
        bed12="{data}.mapped.bam.bed12",
    output:
        named="{data}.mapped.bam.bed12.fasta.named",
    params:
        sd=SD
    shell:"""
{params.sd}/RenameFastaWithGenes.py {input.fa} {input.bed12} > {output.named}
"""

rule MapNamed:
    input:
        fa="{data}.mapped.bam.bed12.fasta.named",
        asm="assembly.orig.fasta"
    output:
        mapped="{data}.mapped.bam.bed12.fasta.named.mm2",
    params:
        grid_opts=config["grid_large"]
    resources:
        load=16
    shell:"""
minimap2 {input.asm} {input.fa} -t 16 > {output.mapped}
"""

rule SelectDups:
    input:
        mm2="{data}.mapped.bam.bed12.fasta.named.mm2",
    output:
        dups="{data}.mapped.bam.bed12.fasta.named.mm2.dups",
    params:
        sd=SD
    shell:"""
{params.sd}/SelectDuplicationsFromMM2.py {input.mm2} > {output.dups}
"""

rule SelectDupsOneIsoform:
    input:
        dups="gencode.mapped.bam.bed12.fasta.named.mm2.dups",
    output:
        iso="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform"
    params:
        sd=SD
    shell:"""
cat {input.dups} | {params.sd}/FilterMembersFromSameIsoformSet.py stdin | {params.sd}/SimplifyName.py | bedtools groupby -g 1,6,8,9 -c 1 -o first -full | cut -f 1-12 | sort > {output.iso}
"""

rule GetGeneCoverage:
    input:
        iso="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform",
        bins="hmm/coverage.bins.bed.gz",
        mean="hmm/mean_cov.txt"
    output:
        cov="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.bed.txt",
        bed="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.bed"
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:"""
cat {input.iso} | awk '{{ print $6"\\t"$8"\\t"$9"\\t"$0;}}' | bedtools groupby -g 1-3 -o first -full -c 1 >  {input.iso}.bed
{params.sd}/GetCoverageOfRegions.sh {input.iso}.bed {input.bins} {input.mean} 
"""

rule GetCombinedTable:
    input:
        bed="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.bed",
        cov="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.bed.txt"
    output:
        combined="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined"
    params:
        sd=SD
    shell:"""
paste {input.bed} {input.cov} | awk '{{ c=int($NF); if (c < 2) {{ c=2; }} print $1"\\t"$2"\\t"$3"\\t"$4"\\t"c-2;}}' > {output.combined}
"""

rule AddCollapsedGenes:
    input:
        comb="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined",
        coll="collapsed_dups_with_genes.bed"
    output:
        colm="collapsed_dups_with_genes.bed.not_multimapped",
        colx="collapsed_dups_with_genes.bed.exclusive",                
    params:
        sd=SD
    shell:"""
{params.sd}/AppendCollapsedDupsList.py {input.comb} {input.coll} > {output.colm}
bedtools intersect -v -a {output.colm} -b {input.comb} | \
   sed "s/collapsed\///g" | \
   awk '{{ $4-=2; if ($4 < 0) {{ $4 = 0; }} print; }}' > {output.colx}
"""

   
rule MinimapGeneModelBed:
    input:
        bam="{data}.mapped.bam"
    output:
        bed="{data}.mapped.bam.bed12"
    params:
        grid_opts=config["grid_large"]
    resources:
        load=1
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
    resources:
        load=1
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
    resources:
        load=1
    shell:"""
head -5 split-rna/{wildcards.data}-1.fa.psl > {output.aln}
for p in {input.psl}; do
tail -n +6 $p 
done  >> {output.aln}
"""


#rule LinkOriginalAssembly:
#    input:
#        ref=assembly
#    output:
#        asm="assembly.orig.fasta",
#        fai="assembly.orig.fasta.fai",
#    shell:"""
#ln -s {input.ref} {output.asm}
#samtools faidx {output.asm}
#"""
    
rule PlotIdeogram:
    input:
        fai="assembly.orig.fasta.fai",
        bed="sedef_out/high_ident/final.sorted.bed.final.filt"
    output:
        pdf=config["species"] + ".sd_dist.pdf"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
Rscript {params.sd}/PlotIdeogram.R {input.fai} {input.bed} {output.pdf}
"""

rule SortDups:
    input:
        comb="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map",
    output:
        comb="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.sorted",
    shell:"""
bedtools sort -i {input.comb} > {output.comb}
"""

rule GetDepthOverDups:
    input:
        comb="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.sorted",
        bins="hmm/coverage.bins.bed.gz",
        avg="hmm/mean_cov.txt"
    output:
        depth="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth",
    shell:"""
m=`cat hmm/mean_cov.txt`

bedtools intersect -loj -a {input.comb} -b {input.bins} -sorted | bedtools groupby -g 1-4 -c 10 -o mean | cut -f 4,5 | awk -v m=$m '{{ print $1"\\t"$2/m;}}' > {input.comb}.cn
paste {input.comb} <( cut -f 2 {input.comb}.cn ) | sort -k4,4 > {output.depth}
"""

rule FilterLowDepthDups:
    input:
        depth="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth",
    output:
        depth_filt="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt",
    shell:"""
cat {input.depth} | awk '{{ if ($NF > 0.05) print;}}' > {output.depth_filt}

"""

rule GetFullGeneCountTable:
    input:
        depth_filt="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt",
    output:
        gene_count="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.gene_count",
        gene_count_2column="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.gene_count_multi_single",
        asm_gene_count="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.asm_gene_count",
    params:
        sd=SD
    shell:"""
cat {input.depth_filt} | awk '{{ $5=1; print;}}' | tr " " "\\t" | bedtools groupby -g 4 -c 5 -o sum | awk '{{ if ($2 > 1) {{ print;}} }}' > {output.asm_gene_count}

{params.sd}/CountGenes.py < {input.depth_filt} > {output.gene_count_2column}

cat {input.depth_filt} | awk '{{ if ($5 == 0) {{ $5=1;}} print;}}' | tr " " "\\t" | bedtools groupby -g 4 -c 5 -o sum | awk '{{ if ($2 > 1) {{ print;}} }}' > {output.gene_count}
"""

rule cramBam:
    input:
        bam=config['bam'],
        orig="assembly.orig.fasta"
    output:
        cram="assembly.cram",
    params:
        grid_opts=config["grid_medium"],
    shell:"""

samtools view {input.bam} -C -@ 4 -T {input.orig} -o {output.cram}
"""

rule RcramBam:
    input:
        rbam="ref_aligned.bam",
        ref="assembly.hg38.fa"
    output:
        rcram="ref_aligned.cram",
    params:
        grid_opts=config["grid_medium"],
    shell:"""
samtools view {input.rbam} -C -@ 4 -T {input.ref} -o {output.rcram}

"""



rule RemoveBams:
    input:
        rbam="ref_aligned.bam",
        rs="hmm_ref/collapsed_duplications.split.bed",
        don="Rhmm.done",
        bam=config['bam'],
        low_cov_tandem_dups="sedef_out/tandem_dups.low_cov.bed",       
        s="collapsed_duplications.split.bed",
        aln=expand("aligned/{b}.bam", b=bamFiles.keys()),
        Raln=expand("ref_aligned/{b}.bam", b=bamFiles.keys()),        
        asm_gene_count="gencode.mapped.bam.bed12.fasta.named.mm2.dups.one_isoform.txt.combined.and_unique_map.depth.filt.asm_gene_count",
    output:
        d="done.done",
    shell:"""
rm {input.aln}
rm {input.Raln}
touch {output}
 
    """
