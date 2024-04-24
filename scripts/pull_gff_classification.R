# Pull isoforms of interest from PB pigeon output (transcript classification and filtering)
# Outputs 1) classification table entries and 2) isoform GFF

setwd("/data/khanlab3/gb_omics/projects/caplen/ccrgb_13/analysis/isoseq_pipeline")

library(plyr)

# Read in gene name, ID, and sample
args = commandArgs(trailingOnly=TRUE)
gene = args[1]
id = args[2]
sample = args[3]

# Load classification table for the sample
classify = read.table(paste("classify/",sample,"_classification.filtered_lite_classification.txt",sep=""),header=TRUE,sep='\t')
print("Loaded classification table")

# Filter to gene of interest
classify.goi = classify[which(classify$associated_gene == id),]
goi.ids = classify.goi$isoform
rm(classify)
print("Filter to gene of interest")

# Load GFF for each sample, filter by PB ID
gff = read.table(paste("classify/",sample,".sorted.filtered_lite.gff",sep=""),sep='\t')
goi = ldply(goi.ids,function(y) gff[grep(paste("transcript_id ",y,";",sep=""),gff$V9),])
rm(gff)

print(paste("Pull GFF"))
  
# Write all
write.table(classify.goi,file=paste("genes_of_interest/",sample,".",gene,"_classification.filtered_lite_classification.txt",sep=""),
            sep='\t',quote=FALSE,row.names=FALSE)
write.table(goi,file=paste("genes_of_interest/",sample,".",gene,".sorted.filtered_lite.gff",sep=""),
            sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
