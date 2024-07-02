# Pull isoforms of interest from PB pigeon output (transcript classification and filtering)
# Outputs 1) classification table entries and 2) isoform GFF

library(plyr)

# Read in gene name, ID, and sample
args = commandArgs(trailingOnly=TRUE)
id = args[1]
classification_file = args[2]
gff_file = args[3]
out_classification = args[4]
out_gff = args[5]

# Load classification table for the sample
classify = read.table(classification_file,header=TRUE,sep='\t')
print("Loaded classification table")

# Filter to gene of interest
classify.goi = classify[which(classify$associated_gene == id),]
print(paste("Entries: ",dim(classify.goi)[1],sep=""))
goi.ids = classify.goi$isoform
rm(classify)
print("Filter to gene of interest")

# Load GFF for each sample, filter by PB ID
gff = read.table(gff_file,sep='\t')
goi = ldply(goi.ids,function(y) gff[grep(paste("transcript_id ",y,";",sep=""),gff$V9),])
rm(gff)
print(paste("Pull GFF"))
  
# Write all
write.table(classify.goi,file=out_classification,sep='\t',quote=FALSE,row.names=FALSE)
write.table(goi,file=out_gff,sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
