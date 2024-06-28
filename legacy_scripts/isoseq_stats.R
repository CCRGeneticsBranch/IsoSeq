# Parse IsoSeq output stats

setwd("~/Genetics Branch/Caplen/Long read quantification")

library(reshape2)
library(openxlsx)
library(ggplot2)

theme_set(theme_bw())

# Load tables
isoseq = read.xlsx("IsoSeq stats 2301614.xlsx",sheet="Isoseq",sep.names=" ")
colnames(isoseq)[1] = "Sample"

pigeon = read.xlsx("IsoSeq stats 2301614.xlsx",sheet="Pigeon",sep.names=" ")
colnames(pigeon)[1] = "Sample"

isoseq = merge(isoseq,pigeon,by="Sample")
rm(pigeon)

isoseq = cbind(isoseq,colsplit(isoseq$Sample,"_",c("X1","X2","siRNA","CellLine")))
isoseq = isoseq[,setdiff(colnames(isoseq),c("X1","X2"))]

# IsoSeq/Pigeon stats
isoseq.reads = melt(isoseq[,c(1:3,5:6,13:14,20:21)],id.vars=c("Sample","siRNA","CellLine"),
                    variable.name="Step",value.name="Count")
isoseq.isoforms = melt(isoseq[,c(1,7,11:12,20:21)],id.vars=c("Sample","siRNA","CellLine"),
                       variable.name="Step",value.name="Count")
isoseq.subcategory = melt(isoseq[,c(1,16:19,20:21)],id.vars=c("Sample","siRNA","CellLine"),
                       variable.name="Category",value.name="Count")
isoseq.genes = melt(isoseq[,c(1,8:10,20:21)],id.vars=c("Sample","siRNA","CellLine"),
                    variable.name="Step",value.name="Count")


png("isoseq_reads.png",width=800)
print(ggplot(isoseq.reads,aes(x=paste(CellLine,siRNA,sep=" "),y=Count,fill=Step)) + geom_bar(stat="identity",position="dodge") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + xlab("Sample") + ylab("Reads"))
dev.off()

ggplot(isoseq.isoforms,aes(x=paste(CellLine,siRNA,sep=" "),y=Count,fill=Step)) + geom_bar(stat="identity",position="dodge") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + xlab("Sample") + ylab("Transcripts")

png("pigeon_category.png",width=800)
print(ggplot(isoseq.subcategory,aes(x=paste(CellLine,siRNA,sep=" "),y=Count,fill=Category)) + geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + xlab("Sample") + ylab("Filtered isoforms"))
dev.off()

ggplot(isoseq.genes,aes(x=paste(CellLine,siRNA,sep=" "),y=Count,fill=Step)) + geom_bar(stat="identity",position="dodge") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + xlab("Sample") + ylab("Genes")

# Saturation
sat = read.xlsx("IsoSeq stats 2301614.xlsx",sheet="Saturation",sep.names=" ")
colnames(sat)[1] = "Sample"
sat = cbind(sat,colsplit(sat$Sample,"_",c("X1","X2","siRNA","CellLine")))
sat = sat[,setdiff(colnames(sat),c("X1","X2"))]

sat = melt(sat,id.vars=c("Sample","siRNA","CellLine","reads"))

png("pacbio_saturation.png",width=800)
print(ggplot(sat,aes(x=reads,y=value,color=paste(CellLine,siRNA,sep=" "))) + geom_line() + facet_wrap(~variable) + 
  xlab("Subsampled reads") + ylab("Count") + scale_color_discrete(name="Sample"))
dev.off()
