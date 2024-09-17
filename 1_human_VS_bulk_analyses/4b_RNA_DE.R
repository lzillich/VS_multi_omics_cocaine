#### Differential expression analysis of RNA-seq data (based on vignette https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) ####
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
seq_info <- read.delim("/path/to/Sample Allocation P2023-033-RNA.txt", header=FALSE)
seq_info <- seq_info[c(1:40),]
seq_info$V2 <- gsub("-CUD","",seq_info$V2)
colnames(seq_info) <- c("rn","Brain_ID")
pheno <- read.csv2("/path/to/pheno/CUD_randomised_20220707_layout.csv")
RIN <- read_excel("/path/to/RNA/DE_results/RNA_info_RIN.xlsx")
RIN <- RIN[c(1:40), c("ID","RIN")]
colnames(RIN) <- c("Brain_ID","RIN")
RIN$RIN <- as.numeric(RIN$RIN)
  
pheno2 <- merge(seq_info,pheno,by="Brain_ID")
pheno2 <- merge(pheno2,RIN,by="Brain_ID")
pheno2$Cocaine_at_death <- 0
pheno2$Cocaine_at_death[grep("ocaine",pheno2$Substance_at_death)] <- 1

# Remove irrelevant pheno data 
pheno3 <- pheno2[,c("Brain_ID","rn","CUD","Age","PMI","pH","RIN","Cocaine_at_death")]
write.table(pheno3,"/path/to/RNA/DE_results/pheno_coc_at_death.txt",sep=";",quote = F,row.names = T,col.names = T)
pheno2 <- pheno2[,c("Brain_ID","rn","CUD","Age","PMI","pH","RIN")]
write.table(pheno2,"/path/to/RNA/DE_results/pheno.txt",sep=";",quote = F,row.names = T,col.names = T)

#mean impute missing pH data 
mean(pheno2$pH,na.rm=T)
pheno2$pH[is.na(pheno2$pH)] <- 6.335
pheno2$CUD[pheno2$CUD==1] <- "CUD"
pheno2$CUD[pheno2$CUD==0] <- "Ctrl"

DT <- data.table
DF <- data.frame

#### 2. Import counts ####

# CUD
mRNA_CUD <- get(load("/path/to/RNA/2_quantified/mRNA_CUD_counts.Rdata"))
colnames(mRNA_CUD$counts) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$counts))
mRNA_CUD$targets <- gsub("Aligned.sortedByCoord.out.bam","",mRNA_CUD$targets)
colnames(mRNA_CUD$stat) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$stat))

mRNA_CUD_counts <- DF(mRNA_CUD$counts)

#### 3. START ANALYSIS ####

#pheno

pheno_mRNA <- pheno2
pheno_mRNA <- pheno_mRNA[order(pheno_mRNA$rn,decreasing = T),]
mRNA_CUD_counts <- mRNA_CUD_counts[,order(colnames(mRNA_CUD_counts),decreasing = T)]

pheno_mRNA$rn  == colnames(mRNA_CUD_counts)
# TRUE, continue 

#coldata
coldata <- as.data.frame(pheno_mRNA)
coldata$CUD <- as.factor(coldata$CUD)
coldata$Age <- scale(coldata$Age)
coldata$PMI <- scale(coldata$PMI)
coldata$RIN <- scale(coldata$RIN)
coldata$pH <- scale(coldata$pH)

rownames(coldata) <- coldata$rn
  
#sort and remove PCA outliers
mRNA_CUD_counts <- mRNA_CUD_counts[,!(colnames(mRNA_CUD_counts) %in% c(	
  "DE42NGSUKBR136478","DE80NGSUKBR136473"))]
coldata <- coldata[colnames(mRNA_CUD_counts),]
identical(colnames(mRNA_CUD_counts),rownames(coldata)) # must be TRUE

## Specify model
dds <- DESeqDataSetFromMatrix(countData = mRNA_CUD_counts,
                              colData = coldata,
                              design = ~ CUD + Age + RIN + pH + PMI)
# Prefiltering
keep <- rowSums(counts(dds) >= 4) >=8
dds <- dds[keep,]

# DE analysis
dds <- DESeq(dds)

## Define contrasts, extract results table 
contrast <-  c("CUD", "CUD", "Ctrl")

res <- results(dds, contrast=contrast, alpha = 0.05)

#Order genes by adjusted p-value
res_genes <- as.data.frame(res)
res_genes_Ordered <- res_genes[order(res_genes$padj,decreasing = F), ]
write.table(res_genes_Ordered,"/path/to/RNA/DE_results/DE_results.txt",quote=F,row.names = T,sep=",")


