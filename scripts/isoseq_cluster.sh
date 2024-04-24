#!/bin/bash

# IsoSeq3 clustering pipeline

SAMPLE=$1
echo $SAMPLE

THREADS=$((SLURM_CPUS_PER_TASK-4))

# Make output directory
if [ ! -d cluster ];then mkdir cluster ;fi

# Lima for primer removal
echo "Lima primer removal"
lima ccs_reads/$SAMPLE\.ccs.bam primers.fasta cluster/$SAMPLE\.fl.bam --isoseq --peek-guess -j $THREADS

# IsoSeq3 refine reads
echo "IsoSeq refine"
isoseq3 refine cluster/$SAMPLE\.fl.NEB_5p--NEB_Clontech_3p.bam primers.fasta cluster/$SAMPLE\.flnc.bam --require-polya -j $THREADS
rm cluster/$SAMPLE\.fl.NEB_5p--NEB_Clontech_3p.bam
rm cluster/$SAMPLE\.fl.NEB_5p--NEB_Clontech_3p.bam.pbi

# IsoSeq3 cluster
echo "IsoSeq cluster"
isoseq3 cluster cluster/$SAMPLE\.flnc.bam cluster/$SAMPLE\.clustered.bam --verbose --use-qvs -j $SLURM_CPUS_PER_TASK
