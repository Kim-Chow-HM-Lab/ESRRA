# Code for Bulk Transcriptomics Analysis
# This set of code was used to analyse all the bulk transcriptomics data in this manuscript. 

setwd("~/DEG")

# Data cleanup
library(dplyr)
library(tidyr)
library(DESeq2)

all <- read.table("RawCounts_Full.txt",header=T, sep="\t")
data <- all %>% select(-c("Chr","Start","End","Strand","Length"))

write.table(data, file="RawCounts.txt",quote=F, sep="\t",row.names=F)

#DESeq2
library(dplyr)
library(tidyr)
library(DESeq2)
library(clusterProfiler)
library(biomaRt)
library(ggplot2)
library(ComplexHeatmap)
library(WGCNA)
library(AnnotationDbi)
library(EnhancedVolcano)
library(org.Mm.eg.db)

data <- read.table("RawCounts.txt",header=T, sep="\t",row.names=1)
meta <- read.table("Meta.txt",header=T, sep="\t",row.names=1)

dds <- DESeqDataSetFromMatrix(countData=data,
                              colData=meta,
                              design=~factor(Strain))

dds
dds <- DESeq(dds)
res <- results(dds)
head(results(dds,tidy=TRUE))
summary <- summary(res)

write.table(res,"DEG_results.txt",quote=F,sep="\t")
normalised_counts <- counts(dds,normalized = TRUE)
write.table(normalised_counts, file="Males_NormalisedCounts.txt", sep="\t", quote=F, col.names=NA)
sum(res$padj<0.05, na.rm=T) #4256

# Volcano Plot
pdf(file = "1_Enhanced Volcano Plot_all.pdf")
EnhancedVolcano(res,
                lab = NA,
                ylab = bquote(~Log[10]~'P-adj'), 
                pCutoff=0.05,
                FCcutoff = 0.25,
                x = 'log2FoldChange',
                y = 'padj', 
                xlim = c(-5,10),
                ylim = c(0,60),
                legendLabels = c('Not sig.', 'Log(base2)FC', 'p-adj','p-adj & Log(base2)FC'))
dev.off()

# PCA Plot
vsdata <- vst(dds, blind=FALSE)
pdf(file="2_PCA Plot_all.pdf")
plotPCA(vsdata,intgroup="Strain")
dev.off()

library(ggplot2)
pdf(file = "2a_PCA Plot Annotated.pdf")
vsdata <- vst(dds,blind=FALSE)
z <- plotPCA(vsdata,intgroup="Strain")
nudge <- position_nudge(y = 1)
z + geom_text(aes(label = name), position = nudge)
dev.off()

