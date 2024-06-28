#!/bin/bash

ml minimap2
ml samtools
ml seqtk
ml R

# Usage: bash get_gene_isoform.sh [Gene name of interest] [Gene ID of interest] [Sample]
# Run from /data/khanlab3/gb_omics/projects/caplen/ccrgb_13/analysis/isoseq_pipeline

# EWSR1 - "ENSG00000182944.17"
# FLI1 - "ENSG00000151702.16"
# CD99 - "ENSG00000002586.19"
# PURPL - "ENSG00000250337.5"

# Load sample name and gene of interest
GENE=$1
ID=$2
SAMPLE=$3

REF="/data/khanlab3/gb_omics/Erica/software/genomes/gencode_human_28/GRCh38.primary_assembly.genome.fa"

# Get pigeon classifications and GFF for gene of interest / sample
Rscript /data/khanlab3/gb_omics/projects/caplen/ccrgb_13/scripts/pull_gff_classification.R $GENE $ID $SAMPLE

# Get full fasta headers per isoform
FASTA="collapse/$SAMPLE.hq.collapsed.fasta"

while read line
	do 
		grep $line"|" $FASTA >> genes_of_interest/$SAMPLE\.$GENE\.tmp
	done < <(cut -f1 genes_of_interest/$SAMPLE\.$GENE\_classification.filtered_lite_classification.txt)
sed -e 's/^>//' genes_of_interest/$SAMPLE\.$GENE\.tmp > genes_of_interest/$SAMPLE\.$GENE\.list
rm genes_of_interest/$SAMPLE\.$GENE\.tmp

# Filter fasta based on isoform ID
seqtk subseq $FASTA genes_of_interest/$SAMPLE\.$GENE\.list > genes_of_interest/$SAMPLE\.$GENE\.fasta
rm genes_of_interest/$SAMPLE\.$GENE\.list
echo "Pull isoforms"

# Map fasta to genome using minimap2
minimap2 -ax splice -t $SLURM_CPUS_PER_TASK -uf --secondary=no -C5 $REF genes_of_interest/$SAMPLE\.$GENE\.fasta | samtools view -b -o genes_of_interest/$SAMPLE\.$GENE\.bam
echo "Map isoforms"

# Sort and index bam
samtools sort -o genes_of_interest/$SAMPLE\.$GENE\.sorted.bam -@ $SLURM_CPUS_PER_TASK genes_of_interest/$SAMPLE\.$GENE\.bam
samtools index genes_of_interest/$SAMPLE\.$GENE\.sorted.bam
rm genes_of_interest/$SAMPLE\.$GENE\.bam
echo "Sort and index"
