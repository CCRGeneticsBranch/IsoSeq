# Author: Erica Pehrsson
# Convert IsoSeq Pigeon output (GTF) to MAJIQ input
# Should match IsoQuant output

library(stringr)

# Read in files
args = commandArgs(trailingOnly=TRUE)
gtf_file = args[1]
classification_file = args[2]
out_gtf = args[3]
out_counts = args[4]

# Read in GTF
gtf = read.table(gtf_file,sep='\t',quote="")
print("Read GTF")

# Read in classification table
classification = read.table(classification_file,sep='\t',header=TRUE)
classification = classification[,c("isoform","structural_category","associated_gene","associated_transcript",
                                   "fl_assoc")]
print("Read classification table")

# Get isoform ID from gtf
gtf$isoform = str_match(gtf$V9, "transcript_id \"(.*?)\";")[,2]

# Add fields to gtf
gtf$associated_gene = classification[match(gtf$isoform,classification$isoform),]$associated_gene
gtf$associated_transcript = classification[match(gtf$isoform,classification$isoform),]$associated_transcript
gtf$V9 = gsub("gene_id","pb_gene_id",gtf$V9)
gtf$attribute = paste(gtf$V9,"; reference_transcript_id \"",gtf$associated_transcript,
                      "\"; gene_id \"",gtf$associated_gene,"\"",sep="")
gtf$attribute = gsub(";;",";",gtf$attribute)
gtf$attribute = ifelse(gtf$V3 == "exon",paste(gtf$attribute,";",sep=""),gtf$attribute)

# Write out gtf
write.table(gtf[,c("V1","V2","V3","V4","V5","V6","V7","V8","attribute")],
            file=out_gtf,
            sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# Write out feature counts
counts = classification[,c("isoform","fl_assoc")]
write.table(counts,file=out_counts,
            sep='\t',row.names=FALSE,quote=FALSE)
