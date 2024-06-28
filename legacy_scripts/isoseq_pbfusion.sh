#!/bin/bash

# PacBio fusion detection

SAMPLE=$1
GTF="gencode.v28.annotation.sorted.gtf"

# Make output directory
if [ ! -d fusion ];then mkdir fusion ;fi

# Detect fusions
pbfusion discover -c 5 -s 10000 -M 0 -t $SLURM_CPUS_PER_TASK --gtf $GTF --output-prefix fusion/$SAMPLE collapse/$SAMPLE\.hq.mapped.bam
