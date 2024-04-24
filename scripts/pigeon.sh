#!/bin/bash

# Pigeon classification and filtering

REF="/data/khanlab3/gb_omics/Erica/software/genomes/gencode_human_28/GRCh38.primary_assembly.genome.fa"
GTF="gencode.v28.annotation.sorted.gtf"

SAMPLE=$1
echo $SAMPLE

# Make output directory
if [ ! -d classify ];then mkdir classify ;fi

# Sort transcript GFF
echo "Sort GFF"
pigeon sort collapse/$SAMPLE\.hq.collapsed.gff -o classify/$SAMPLE\.sorted.gff

# Classify isoforms
echo "Classify isoforms"
pigeon classify classify/$SAMPLE\.sorted.gff $GTF $REF -j $SLURM_CPUS_PER_TASK -d "classify" -o $SAMPLE \
    --gene-id --fl collapse/$SAMPLE\.hq.collapsed.abundance.txt

# Filter isoforms
echo "Filter isoforms"
pigeon filter classify/$SAMPLE\_classification.txt --isoforms classify/$SAMPLE\.sorted.gff -j $SLURM_CPUS_PER_TASK

# Report gene saturation
echo "Report gene saturation"
pigeon report classify/$SAMPLE\_classification.filtered_lite_classification.txt classify/$SAMPLE\.saturation.txt -j $SLURM_CPUS_PER_TASK
