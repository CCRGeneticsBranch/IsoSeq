#!/bin/bash

# Sort and index genome annotation file
GTF=/data/khanlab3/gb_omics/Erica/software/gencode/gencode.v28.annotation.gtf

pigeon sort $GTF -o gencode.v28.annotation.sorted.gtf
pigeon index gencode.v28.annotation.sorted.gtf
