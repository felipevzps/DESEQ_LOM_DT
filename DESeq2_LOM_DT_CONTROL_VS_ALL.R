if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db")

install.packages("jsonlite")
install.packages("WGCNA")
install.packages("pheatmap")

library(DESeq2)
library(tximport)
library(jsonlite)
library(pheatmap)

# Reset environment variables
rm(list=ls())

setwd("/home/felipevzps/Documentos/DESEQ_LOM_DT")

targets<-read.csv("target_without_3DT_1Control_LOM.csv",header=TRUE)
rownames(targets)<-targets$SampleName
targets

files <- paste("/home/felipevzps/Documentos/DESEQ_LOM_DT/quant/", targets$SampleName, "/quant.sf",sep="")
tx2gene<-read.delim("transcript2gene_atualizado.txt",header=FALSE)
names(files) <- targets$SampleName
all(file.exists(files))

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
head(txi.salmon$abundance)

Chla_TPM<-txi.salmon$abundance
txi.salmon$abundance

dds <- DESeqDataSetFromTximport(txi.salmon, colData = targets, design = ~Condition)
dds$Condition <-relevel(dds$Condition, ref='CONTROL')

#Merge Technical Replicates
ddsColl <- collapseReplicates(dds, dds$Sample, dds$Run)
colData(ddsColl)
colnames(ddsColl)
ddsColl$runsCollapsed

head(ddsColl)
ddsColl <- DESeq(ddsColl)

vsd <- vst(ddsColl, blind=FALSE)
vsd
####

#dds <- DESeq(dds)

#Check distribution of samples in a PCA, showing the top 500 varying genes
plotPCA(vsd, intgroup=c("Time"),ntop=500)
plotPCA(vsd, intgroup=c("Condition"),ntop=500)
plotPCA(vsd, intgroup=c("Condition","Time"),ntop=500)

### including name to samples in plot
install.packages("ggrepel")

library(ggplot2)
library(ggrepel)

plotPCA(vsd, intgroup=c("Condition","Time"),ntop=500) + geom_text(aes(label=name), vjust=2) + geom_point(size=1)

#dev.copy(png, "PCA_LOM_DT.png")
#dev.off()

###

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Time, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Time, sep="-")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

res_SALT_vs_CONTROL <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs", contrast = c('Condition','CONTROL','SALT'), alpha=0.05)
summary(res_SALT_vs_CONTROL)

sig_SALT_vs_CONTROL<-res_SALT_vs_CONTROL[which(res_SALT_vs_CONTROL$padj<0.05),]
summary(sig_SALT_vs_CONTROL)
dim(sig_SALT_vs_CONTROL)
head(sig_SALT_vs_CONTROL,20)

#Exploring DGE results
drawLines <- function() abline(h=c(-1,1),col="red",lwd=2)
plotMA(res_SALT_vs_CONTROL); drawLines()

TPM_SALT_vs_CONTROL<-Chla_TPM[which(rownames(Chla_TPM) %in% rownames(sig_SALT_vs_CONTROL)),
                              as.character(targets[which(targets$Condition %in% c('SALT','CONTROL')),'SampleName'])]

annot_col= data.frame(Time=factor(targets$Time),
                      Condition=factor(targets$Condition))
annot_col
rownames(annot_col)<-targets$SampleName
pheatmap(TPM_SALT_vs_CONTROL, scale='row', annotation_col = annot_col)

