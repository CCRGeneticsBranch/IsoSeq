#!/bin/bash

# IsoSeq3 mapping and collapse

SAMPLE=$1
REF="/data/khanlab3/gb_omics/Erica/software/genomes/gencode_human_28/GRCh38.primary_assembly.genome.fa"

# Make output directory
if [ ! -d collapse ];then mkdir collapse ;fi

# Map transcripts to reference genome
echo "Map transcripts"
pbmm2 align --preset ISOSEQ --sort cluster/$SAMPLE\.clustered.hq.bam $REF collapse/$SAMPLE\.hq.mapped.bam -j $SLURM_CPUS_PER_TASK

# Collapse redundant transcripts into unique isoforms via exon structure
echo "Collapse transcripts"
isoseq3 collapse --do-not-collapse-extra-5exons collapse/$SAMPLE\.hq.mapped.bam collapse/$SAMPLE\.hq.collapsed.gff -j $SLURM_CPUS_PER_TASK
