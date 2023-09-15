import os
import tempfile
import subprocess
import os.path


# Config
configfile: "sd_analysis.json"

if "asm" not in config and "assembly" in config:
    config["asm"] = config["assembly"]

assembly="assembly.orig.fasta"
asmFai="assembly.orig.fasta.fai"
#asmFai=config["asm"] + ".fai"
asmFaiFile=open(asmFai)
contigs=[l.split()[0] for l in asmFaiFile]
strIdx=[str(i) for i in range(0,len(contigs))]

tempDir=config["temp"]
# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

rule all:
    input:
        mask="assembly.repeat_masked.fasta",
        maskedGenomeOut="assembly.repeat_masked.fasta.out"        


rule SplitContig:
    input:
        asm=assembly
    output:
        contig="split/to_mask.{index}.fasta",
    params:
        contigName=lambda wildcards: contigs[int(wildcards.index)],
        grid_opts=config["grid_small"]
    shell:"""
mkdir -p split
samtools faidx {input.asm} \"{params.contigName}\" > {output.contig}
"""

rule MaskContig:
    input:
        contig="split/to_mask.{index}.fasta"
    output:
        mask="masked/to_mask.{index}.fasta.masked",
        maskOut="masked/to_mask.{index}.fasta.out"
    params:
        grid_opts=config["grid_repeatmasker"],
        repeatLibrary=config["repeat_library"],
        sd=SD        
    shell:"""
mkdir -p masked
TEMP="$TMPDIR/$$_$RANDOM/"
mkdir -p $TEMP
cp \"{input.contig}\" \"$TEMP/to_mask.{wildcards.index}.fasta\" && \
pushd $TEMP &&  \
RepeatMasker {params.repeatLibrary} -pa 8  -s -xsmall \"to_mask.{wildcards.index}.fasta\" && \
popd && \
if [ ! -e $TEMP/to_mask.\"{wildcards.index}\".fasta.masked ]; then
  cp split/to_mask.\"{wildcards.index}\".fasta masked/to_mask.\"{wildcards.index}\".fasta.masked
  cp {params.sd}/repeat_masker.out.header masked/to_mask.\"{wildcards.index}\".fasta.out
else
  cp $TEMP/to_mask.\"{wildcards.index}\".fasta.* masked/ || true
fi
#rm -rf $TEMP
"""

rule SpecialMaskContig:
    input:
        mask="masked/to_mask.{index}.fasta.masked",
    output:
        tt="t2t/to_mask.{index}.fasta.masked",
        ttOut="t2t/to_mask.{index}.fasta.out",
    params:
        grid_opts=config["grid_repeatmasker"],
        repeatLibrary=config["t2t_repeat_library"],
        sd=SD
    shell:"""
if [ {params.repeatLibrary} != "na" ]
then
  mkdir -p t2t
  TEMP="$TMPDIR/$$_$RANDOM/"
  mkdir -p $TEMP
  # 
  # Copy input file to temp dir that should have fast IO
  #
  {params.sd}/hardmask {input.mask} $TEMP/to_mask.{wildcards.index}.fasta && \
  pushd $TEMP &&  \
  RepeatMasker -nolow -libdir /home1/mchaisso/miniconda3/share/RepeatMasker/Libraries/Appended -species human -pa 8 -s -xsmall \"to_mask.{wildcards.index}.fasta\" && \
  popd && \

  if [ ! -e $TEMP/to_mask.\"{wildcards.index}\".fasta.masked ]; then
    ls -l masked/to_mask.\"{wildcards.index}\".fasta.masked
    cp masked/to_mask.\"{wildcards.index}\".fasta.masked t2t/to_mask.\"{wildcards.index}\".fasta.masked
    head -3 masked/to_mask.\"{wildcards.index}\".fasta.out > t2t/to_mask.\"{wildcards.index}\".fasta.out
  else
    cp $TEMP/to_mask.\"{wildcards.index}\".fasta.* t2t/ || true
  fi
  rm -rf $TEMP
else
  echo "*** rule SpecialMaskContig made no change (skipped) ***" > t2t/to_mask.\"{wildcards.index}\".fasta.out
  cp {input.mask} t2t/to_mask.\"{wildcards.index}\".fasta.masked >> t2t/to_mask.\"{wildcards.index}\".fasta.out || true
fi
"""

rule MergeMaskerRuns:
    input:
        humLib="masked/to_mask.{index}.fasta.masked",
        humLibOut="masked/to_mask.{index}.fasta.out",
        t2tLib="t2t/to_mask.{index}.fasta.masked",
        t2tLibOut="t2t/to_mask.{index}.fasta.out"        
    output:
        comb="comb/to_mask.{index}.fasta.masked",
        combOut="comb/to_mask.{index}.fasta.out",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:"""
{params.sd}/comask {output.comb} {input.humLib} {input.t2tLib}
echo {input.humLibOut} > comb/to_mask.{wildcards.index}.names
echo {input.t2tLibOut} >> comb/to_mask.{wildcards.index}.names    
{params.sd}/RepeatMasking/AppendOutFile.py {output.combOut} comb/to_mask.{wildcards.index}.names
"""

rule CombineMask:
    input:
        maskedContigs=expand("comb/to_mask.{index}.fasta.masked", index=strIdx),
        maskedContigsOut=expand("comb/to_mask.{index}.fasta.out", index=strIdx)        
    output:
        maskedGenome="assembly.repeat_masked.fasta",
        maskedGenomeOut="assembly.repeat_masked.fasta.out"        
    params:
        grid_opts=config["grid_small"],
        outFileNames=" ".join(["comb/to_mask.{}.fasta.out".format(i) for i in strIdx]),
        inFastaNames=" ".join(["comb/to_mask.{}.fasta.masked".format(i) for i in strIdx]),
        sd=SD
    shell:"""
cat {params.inFastaNames} > {output.maskedGenome}
{params.sd}/RepeatMasking/AppendOutFile.py {output.maskedGenomeOut} {params.outFileNames}
"""
        
