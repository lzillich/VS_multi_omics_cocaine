## Data Preparation MOFA - adapted from Tutorial code https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html and https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html
# Authors: Lea Zillich and Eric Zillich 
# last change: EZ 2024-09-12

library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(data.table)
library(ggplot2)
library(tidyverse)
library(psych)
library(ggpubr)
library(limma)
library(biomaRt)

DF <- data.frame

dir0 <- "/path/to/protein/data_analysis_results_V1/MOFA"
setwd(dir0)

# pheno file 
pheno_mRNA <- read.csv("/path/to/RNA/DE_results/pheno.txt", sep=";")
colnames(pheno_mRNA)[2] <- "rn_mRNA"

pheno_miRNA <- read.csv("/path/to/miRNA/DE_results/pheno.txt", sep=";")
pheno_miRNA <-pheno_miRNA[,c(1,2)]
colnames(pheno_miRNA)[2] <- "rn_miRNA"

pheno_all <- merge(pheno_mRNA,pheno_miRNA,by="Brain_ID")
pheno_all$pH[is.na(pheno_all$pH)] <- 6.335

# Import protein pheno file 
pheno_prot <- read.csv("/path/to/protein/data_analysis_results_V1/pheno_merged_with_proteomics.txt", sep=";")
colnames(pheno_prot)[1] <- "Brain_ID"
pheno_all_prot <- merge(pheno_all, pheno_prot,by="Brain_ID",all.x=T,all.y=T)
pheno_all_prot2 <-pheno_all_prot[,c(1:11,17)]

write.table(pheno_all_prot2,"/path/to/protein/data_analysis_results_V1/MOFA/pheno_all_prot.txt",sep=";",row.names = F,quote=F)

##### miRNA data #####
norm_counts_miR = get(load("/path/to/miRNA/DE_downstream/WGCNA/input_VS_expression_filtered.Rdata"))

#scale and center miRNAseq data
norm_counts_miR <- jyluMisc::mscale(norm_counts_miR, center = TRUE, scale = TRUE, useMad = FALSE)
  colnames(norm_counts_miR) <- pheno_all_prot$Brain_ID[match(colnames(norm_counts_miR),pheno_all_prot$rn_miRNA)]

# Add NA for unavailable samples
pheno_all_prot2$Brain_ID[!(pheno_all_prot2$Brain_ID %in% colnames(norm_counts_miR))]

mat<- matrix(nrow=1542,ncol=3,NA) 
colnames(mat)<- c("ID1","ID2","ID3")
norm_counts_miR_all <- cbind(norm_counts_miR,mat)

####### mRNA data #######
  norm_counts_mRNA = get(load("/path/to/RNA/DE_downstream/WGCNA/input_VS_expression_filtered.Rdata"))
  
  #filter out low variable genes
  ntop=4270
  sds <- genefilter::rowSds(norm_counts_mRNA)
  exprMat<-norm_counts_mRNA[order(sds,decreasing = T)[1:ntop],]
  
  #scale and center RNAseq data
  exprMat <- jyluMisc::mscale(exprMat, center = TRUE, scale = TRUE, useMad = FALSE)
  colnames(exprMat) <- pheno_all_prot$Brain_ID[match(colnames(exprMat),pheno_all_prot$rn_mRNA)]
  
  # Add NA for unavailable samples
  pheno_all_prot2$Brain_ID[!(pheno_all_prot2$Brain_ID %in% colnames(exprMat))]
  mat<- matrix(nrow=4270,ncol=3,NA) 
  colnames(mat)<- c("ID1","ID2","ID3")
  exprMat_all <- cbind(exprMat,mat)
  
####### protein data #######
  norm_prot = read.csv("/path/to/protein/data_analysis_results_V1/normalized_data.txt", row.names=1, sep=";")

  #scale and center protein data
  norm_prot <- jyluMisc::mscale(norm_prot, center = TRUE, scale = TRUE, useMad = FALSE)
  colnames(norm_prot) <- pheno_all_prot2$Brain_ID[match(colnames(norm_prot),pheno_all_prot2$ID)]
  
  # Add NA for unavailable samples
  pheno_all_prot2$Brain_ID[!(pheno_all_prot2$Brain_ID %in% colnames(norm_prot))]
 mat<- matrix(nrow=4270,ncol=1,NA) 
 colnames(mat)<- c("ID4")
 norm_prot_all <- cbind(norm_prot,mat)
 
 list_mofa <- list( norm_counts_miR_all, exprMat_all,norm_prot_all)
  
  save(list_mofa, file ="/path/to/protein/data_analysis_results_V1/MOFA/VS_mRNA_miRNA_protein.Rdata")
  

