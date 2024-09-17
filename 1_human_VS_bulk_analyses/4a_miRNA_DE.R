#### Differential expression analysis of miRNA-seq data (based on DESeq2 vignette https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) ####
# Author: Eric Zillich
# Last change 2024-09-12 EZ

library(ggplot2)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(data.table)
library(readxl)
library(readr)
library(varhandle)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(qqman)


#### 1. Import phenotype data ####
seq_info <- read.delim("/path/to/Sample Allocation P2023-034-miRNA.txt", header=FALSE)
colnames(seq_info) <- c("rn","Brain_ID")
pheno <- read.csv2("/path/to/pheno/CUD_randomised_20220707_layout.csv")
RIN <- read_excel("/path/to/miRNA/DE_results/RNA_info_RIN.xlsx")
RIN <- RIN[c(1:40), c("ID","RIN")]
colnames(RIN) <- c("Brain_ID","RIN")
RIN$RIN <- as.numeric(RIN$RIN)
  
pheno2 <- merge(seq_info,pheno,by="Brain_ID")
pheno2 <- merge(pheno2,RIN,by="Brain_ID")

# Remove irrelevant pheno data 
pheno2 <- pheno2[,c("Brain_ID","rn","CUD","Age","PMI","pH","RIN")]
write.table(pheno2,"/path/to/miRNA/DE_results/pheno.txt",sep=";",quote = F,row.names = T,col.names = T)

#mean impute missing pH data 
mean(pheno2$pH,na.rm=T)
pheno2$pH[is.na(pheno2$pH)] <- 6.335

DT <- data.table
DF <- data.frame

#### 2. Import counts ####

# CUD
miRNA_CUD <- get(load("/path/to/miRNA/2_quantified/miRNA_CUD_counts.Rdata"))
colnames(miRNA_CUD$counts) <- gsub("_1_trimmedAligned.sortedByCoord.out.bam","",colnames(miRNA_CUD$counts))
miRNA_CUD$targets <- gsub("_1_trimmedAligned.sortedByCoord.out.bam","",miRNA_CUD$targets)
colnames(miRNA_CUD$stat) <- gsub("_1_trimmedAligned.sortedByCoord.out.bam","",colnames(miRNA_CUD$stat))

miRNA_CUD_counts <- DF(miRNA_CUD$counts)


#### 3. START ANALYSIS ####

#pheno

pheno_miRNA <- pheno2
pheno_miRNA <- pheno_miRNA[order(pheno_miRNA$rn,decreasing = T),]


miRNA_CUD_counts <- miRNA_CUD_counts[,order(colnames(miRNA_CUD_counts),decreasing = T)]
miRNA_CUD_counts <- miRNA_CUD_counts

pheno_miRNA$rn  == colnames(miRNA_CUD_counts)
# TRUE, continue 

#coldata
coldata <- as.data.frame(pheno_miRNA)
coldata$CUD <- as.factor(coldata$CUD)
coldata$Age <- scale(coldata$Age)
coldata$PMI <- scale(coldata$PMI)
coldata$RIN <- scale(coldata$RIN)
coldata$pH <- scale(coldata$pH)

rownames(coldata) <- coldata$rn

#sort and remove PCA outliers
miRNA_CUD_counts <- miRNA_CUD_counts[,!(colnames(miRNA_CUD_counts) %in% c(	
  "DE42NGSUKBR136478","DE80NGSUKBR136473"))]
coldata <- coldata[colnames(miRNA_CUD_counts),]
identical(colnames(miRNA_CUD_counts),rownames(coldata)) # must be TRUE

## Specify model
dds <- DESeqDataSetFromMatrix(countData = miRNA_CUD_counts,
                              colData = coldata,
                              design = ~ CUD + Age + RIN + pH + PMI)
# Prefiltering
keep <- rowSums(counts(dds) >= 2) >= 4
dds <- dds[keep,]

# DE analysis
dds <- DESeq(dds)

## Define contrasts, extract results table 
contrast <-  c("CUD", "1", "0")

res <- results(dds, contrast=contrast, alpha = 0.05)
res_genes <- as.data.frame(res) 
res_genes_Ordered <- res_genes[order(res_genes$padj,decreasing = F), ]
write.table(res_genes_Ordered,"/path/to/miRNA/DE_results/DE_results.txt",quote=F,row.names = T,sep=",")

