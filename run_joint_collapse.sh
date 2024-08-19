#!/bin/bash

pbmm2 align --preset ISOSEQ --sort ../isoseq_updated/results/clustered.bam /data/khanlab3/gb_omics/resources/genomes/gencode_human_28/GRCh38.primary_assembly.genome.fa ../isoseq_updated/results/clustered.mapped.bam -j 24
isoseq3 collapse --do-not-collapse-extra-5exons ../isoseq_updated/results/clustered.mapped.bam ../isoseq_updated/results/clustered.mapped.gff -j 24
