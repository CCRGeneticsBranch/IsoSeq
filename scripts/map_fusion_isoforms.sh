#!/bin/bash

# Map fusion isoforms to custom genome
# Upstream: run PBFusion on IsoSeq clustered transcripts (HQ), get list of EWSR1 fusion transcript IDs
# Input: sample name, custom fusion genome (fasta)
# Output: fusion isoforms mapped to fusion genome (sorted, indexed bam)

ml minimap2
ml samtools
ml seqtk

# Usage: bash map_fusion_isoforms.sh [sample] [Custom fusion reference]

# Get variables
SAMPLE=$1
REF=$2 #../fusion_genome/ewsr1_fli1.fa

# IsoSeq clustered HQ fasta
FASTA=cluster/$SAMPLE\.clustered.hq.fasta.gz

# Fusion isoform IDs
IDS=fusion/bam/$SAMPLE\_ewsr1.list

# Fusion genome
GENOME=$(basename $REF .fa)

# Filter fasta file based on isoform ID
seqtk subseq $FASTA $IDS > fusion/bam/$SAMPLE\.fasta
echo "Pull isoforms"

# Map fasta to genome using minimap2
minimap2 -ax splice -t $SLURM_CPUS_PER_TASK -uf --secondary=no -C5 $REF fusion/bam/$SAMPLE\.fasta | samtools view -b -o fusion/bam/$SAMPLE\_$GENOME\.bam
echo "Map isoforms"

# Sort and index
samtools sort -o fusion/bam/$SAMPLE\_$GENOME\.sorted.bam fusion/bam/$SAMPLE\_$GENOME\.bam
samtools index fusion/bam/$SAMPLE\_$GENOME\.sorted.bam
rm fusion/bam/$SAMPLE\_$GENOME\.bam
echo "Sort and index"
