import os
import tempfile
import subprocess
import os.path


# Config
configfile: "sd_analysis.json"




tempp=config['temp']
if "temp2" not in config:
    config["temp2"] = config["temp"]
    
if config['temp2']!="":
    tempp=config['temp2']


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

assembly="assembly.orig.fasta"


subs=["all", "high_ident"]


localrules: all, AnnotateResolvedTandemDups, GetUniqueGencodeUnresolvedDupGenes,  IntersectGenesWithFullSDList, FullDupToBed12, FullDupToLinks, MakeWMBed, MaskFile, ConvertHMMCopyNumberToCollapsedDuplications, SortSedef, FilterSedef, CountMaskedSedef, RemoveSedefTooMasked, MakeSedefGraph, MakeSedefGraphTable, FilterByGraphClusters, FullDupToBed12, FiltDupToBed12, GetUniqueGencodeUnresolvedDupGenesCN, GetUniqueGencodeUnresolvedDupGenes, GetGencodeMulticopy, GetGencodeMappedInDup, GetSupportedMulticopy,FindResolvedDuplicatedGenes, Bed12ToBed6, CombineGenesWithCollapsedDups, CombineDuplicatedGenes, MinimapGeneModelBed, FilterGencodeBed12, FindGenesInResolvedDups, SelectOneIsoform, SplitSplicedAndSingleExon, AnnotateLowCoverageFlanks, UnionMasked,GetNamedFasta, SelectDups, SortDups, GetDepthOverDups, FilterLowDepthDups, GetFullGeneCountTable, AddCollapsedGenes, GetCombinedTable, SelectDupsOneIsoform, GetFinalMerged, DupsPerContig, GetAllMultiGenes, AnnotateHighIdentity, GetTotalMasked, AnnotateResolvedTandemDups, GeneCountFact, GetFullGeneCountTable, FilterMultiExonBed, MappedSamIdentityDups, RemoveOriginal, RemoveBams, MakeSedefIntv, HighestIdentPairs, SelectHighIdent, GetCollapseByRange, GetCollapsedMask, GetCN



config["assembly"] = "assembly.repeat_masked.fasta"

#import shutil
#onsuccess:
#    shutil.rmtree(".snakemake")
subs=["all", "high_ident"]

if "name_map" not in config:
    config["name_map"] = "NO_OP"

rule all:
    input:
        fai=assembly+".fai",
        sedef="sedef_out/final.bed",
        sedef_sorted="sedef_out/final.sorted.bed",
        filt="sedef_out/all/final.sorted.bed.final.filt",        
        asmMask=expand("{asm}.count_masked", asm=["assembly.orig.fasta", "assembly.masked.fasta", "assembly.union_masked.fasta"]),



#
# Simple preprocessing, make sure there is an index on the assembly.
#


rule MakeFaiLinkOrig:
    input:
        asm=config['assembly']
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

def GetBam(f):
    return f.split("/")[-1]

#
# Map individual bams separately
#

#rule MakeWMDB:
#    input:
#        asm=assembly
#    output:
#        wm_db="wmdb"
#    params:
#        grid_opts=config["grid_medium"]
#    resources:
#        load=2
#    shell:"""
#windowmasker -mk_counts -in {input.asm} -out {output.wm_db}  || true
#"""
#
#rule MakeWMIntv:
#    input:
#        wm_db="wmdb",
#        asm=assembly
#    output:
#        intv="wm_mask_intervals"
#    params:
#        grid_opts=config["grid_large"]
#    resources:
#        load=2
#    shell:"""
#windowmasker -ustat {input.wm_db} -in {input.asm} -out {output.intv}  || true
#"""
#
#rule MakeWMBed:
#    input:
#        intv="wm_mask_intervals"
#    output:
#        bed="wm_mask_intervals.bed"
#    params:
#        grid_opts=config["grid_small"]
#    resources:
#        load=1
#    shell:"""
#cat {input.intv} | awk '{{ if (substr($1,0,1) == ">") {{ name=substr($1,1); }} else {{ if ($3-$1 > 100) print name"\\t"$1"\\t"$3;}} }}' | tr -d ">" > {output.bed}
#"""
#
#rule MaskFile:
#    input:
#        bed="wm_mask_intervals.bed",
#        asm=assembly
#    output:
#        masked="assembly.masked.fasta"
#    params:
#        grid_opts=config["grid_medium"],
#        sd=SD
#    resources:
#        load=1
#    shell:"""
#{params.sd}/bemask {input.asm} {input.bed} {output.masked}
#"""
#
#
# Run repeat masker on the assembly. This will be combined with the
# windowmasker to generate a masked genome.
#

    

#
#  The final masked genome combines wm and repeatmasker
#

#rule UnionMasked:
#    input:
#        orig="assembly.orig.fasta",
#        wm="assembly.masked.fasta",
#    output:
#        comb="assembly.union_masked.fasta"
#    params:
#        sd=SD,
#        grid_opts=config["grid_large"]
#    resources:
#        load=1
#    shell:"""
#    {params.sd}/comask {output.comb} {input.orig} {input.wm}
#    samtools faidx {output.comb}
#"""

#
# The following rules do the initial resolved repeat detection with
# sedef, and then postprocess the output to remove excess duplications.
# 

#
# Initial run of sedef.
#
rule RunSedef:
    input:
#        asm="assembly.union_masked.fasta"
    output:
        done="sedef_out/final.bed"
    params:
        asm="assembly.union_masked.fasta",        
        grid_opts=config["grid_sedef"],
        sd=SD
    resources:
        load=8
    shell:"""

module load gcc/8.3.0
module load time
module load parallel
export PATH=$PATH:{params.sd}/sedef

{params.sd}/sedef.sh  {params.asm} -j 8
"""
#
# Sort by chrom and start, fixing a bug in sedef output that misses
# column output.
#

rule SortSedef:
    input:
        bed="sedef_out/final.bed"
    output:
        s="sedef_out/final.sorted.bed",
        done="sedef.done",
    params:
        grid_opts=config["grid_medium"]
    resources:
        load=1
    shell:"""
first=`head -1 {input.bed} | awk '{{ print NF;}}'`
sort -k1,1 -k2,2n {input.bed} | \
  awk -vf=$first '{{ if (NF == f) print;}}' | \
  bedtools groupby -g 1-6 -o first -full -c 1 | cut -f 1-$first > {output.s}

touch {output.done}
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



rule GetCN:
    input:
        vcf="hmm/copy_number.vcf"
    output:
        bed="hmm/copy_number.bed.gz",
    params:
        sd=SD        
    shell:"""
{params.sd}/CovVcfToBed.py {input.vcf} | gzip -c > {output.bed}
"""
    

rule GencodeCN:
    input:
        gc="gencode.mapped.bam.bed12",
        cn="hmm/cov_bins.bed.gz",
        genome="assembly.orig.fasta"
    params:
        grid_opts=config["grid_small"],
    output:
        gccn="gencode.mapped.bam.bed12.cn",
    shell:"""
bedtools intersect -loj -g {input.genome}.fai -a gencode.mapped.bam.bed12 -b {input.cn} -sorted | awk '{{ if (NF== 17) print;}}' | bedtools groupby -g 1-4 -c 17 -o mean > {output.gccn}

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




rule GeneDupsToLinks:
    input:
        bed="gencode.mapped.bam.bed12.multi_exon.fasta.named.mm2.dups",
    output:
        links="circos/genes_in_resolved_dups.links.tsv",
        names="circos/genes_in_resolved_dups.links.names.tsv",
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    resources:
        load=1
    shell:"""
mkdir -p circos

cat  {input.bed}  | \
 grep -v "Gm" | \
 {params.sd}/FilterMembersFromSameIsoformSet.py  stdin | \
 tr "/" "\\t" | \
 {params.sd}/SimplifyName.py  | \
 tr "\:\-" "\\t\\t" | \
 {params.sd}/DupsToArcs.py {output.links} {output.names}

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

rule MakeCircOS:
    input:
        asm="assembly.orig.fasta",
        coll="gencode.mapped.bam.bed12.dups.unique",
        links="circos/genes_in_resolved_dups.links.tsv",
        names="circos/genes_in_resolved_dups.links.names.tsv"
    output:
        plt="circos/circos.png"
    params:
        sd=SD,
        grid_opts=config["grid_medium"],
        name_map=config["name_map"]
    resources:
       load=1
    shell:"""
mkdir -p circos
cat {input.coll} | grep -v "Gm" | tr "#" "_" | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$(NF-1)"\\t"$NF;}}' > {input.coll}.name.cn
cat {input.asm}.fai | tr "#" "_" > {input.asm}.fai.u
cut -f 1 {input.asm}.fai.u | sed "s/_/__/g" | tr "#" "_" > {input.asm}.display_name
paste {input.asm}.fai.u {input.asm}.display_name | awk '{{ if ($2 > 500000) {{ print "chr\\t-\\tvar"$1"\\t"$6"\\t"0"\\t"$2"\\t"$1;}}  }}' | \
 {params.sd}/ReplaceChromName.py NO_OP >  circos/karyotype.txt

cat {input.coll}.name.cn | awk '{{ print "var"$1"\\t"$2"\\t"$3"\\t"$4;}}' | {params.sd}/ReplaceChromName.py {params.name_map} > circos/cn.lab.txt
cat {input.coll}.name.cn | awk '{{ print "var"$1"\\t"$2"\\t"$3"\\t"$5;}}' | {params.sd}/ReplaceChromName.py {params.name_map} > circos/cn.txt

#{params.sd}/MakeDup.py --bed {input.coll}  --collapsed {input.coll}.name.cn --links circos/resolved_dups.txt --labels circos/resolved_dups.labels.txt
cd circos && circos --conf {params.sd}/circos.conf 
"""


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
        dups="sedef_out/all/final.sorted.bed.final.filt"
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
rule FindResolvedDuplicatedGenes:
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
        dups="collapsed_duplications.bed.range4",
        asm="assembly.union_masked.fasta"
    output:
        rnabedout="{data}.mapped.bam.bed12.dups",
    params:
        grid_opts=config["grid_small"],
        sd=SD, 
    resources:
        load=1
    shell:"""


#
# The following intersects mapped genes with collapsed duplications. It also makes sure
# that the number of fields in the intersection add up to the correct toal
# since missing values were causing problems in downstream analysis. 
#

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


#rule CountMaskedAsmp:
#    input:
#        asm="{assembly}"
#    output:
#        countMasked="{assembly}.count_masked"
#    params:
#        grid_opts=config["grid_small"],
#        sd=SD
#    resources:
#        load=1
#    shell:"""
#cat {input.asm} | {params.sd}/nl > {output.countMasked}
#"""
#




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

rule MakeSedefIntv:
    input:
        sedefFinal="sedef_out/all/final.sorted.bed.ident.pairs",
    output:
        sedefIntv="sedef_out/all/final.sorted.bed.intv",
    params:
        sd=SD
    shell:"""
{params.sd}/DivideBedToNonoverlappingIntervals.py {input.sedefFinal} | bedtools sort > {output.sedefIntv}
"""

rule HighestIdentPairs:
    input:
        sedefIntv="sedef_out/all/final.sorted.bed.final.filt",
    output:
        sedefPairs="sedef_out/all/final.sorted.bed.ident.pairs"
    shell:"""
cat {input.sedefIntv} | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$21; print $4"\\t"$5"\\t"$6"\\t"$21;}}' | bedtools sort > {output.sedefPairs}
"""

rule SelectHighIdent:
    input:
        sedefIntv="sedef_out/all/final.sorted.bed.intv",
        sedefPairs="sedef_out/all/final.sorted.bed.ident.pairs"
    output:    
        sedef_high_uniq="sedef_out/all/final.sorted.bed.uniq.high",
    shell:"""
bedtools intersect -a {input.sedefIntv} -b {input.sedefPairs} -loj -sorted | bedtools groupby -g 1-3 -c 7 -o max > {output.sedef_high_uniq}
"""    
