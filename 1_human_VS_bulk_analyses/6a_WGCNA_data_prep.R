## Data Preparation WGCNA 
## Code based on WGCNA tutorials by Langfelder & Horvath
# Adapted by Lea Zillich and Eric Zillich 
# last change: EZ 2024-09-12

library(WGCNA)
library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(DESeq2)

DF <- data.frame

#### miRNA-seq data ####

dir0 <- "path/to/miRNA"

setwd(dir0)
options(stringsAsFactors = FALSE)

#read pheno
pheno_c<- read.csv("/path/to/RNA/DE_results/pheno_coc_at_death.txt", sep=";")

phe_orig <- read.csv2("/path/to/pheno/CUD_randomised_20220707_layout.csv")

#Subset for relevant pheno information
phe_orig <-phe_orig[,c(1,13,14)]
pheno_all <- merge(pheno_c,phe_orig,by="Brain_ID")
pheno_all$pH[is.na(pheno_all$pH)] <- 6.335

# Remove PCA outliers
pheno_all <- pheno_all[!(pheno_all$Brain_ID %in% c("ID1","ID2")),]
pheno_cc <- pheno_all[,-c(9,10)]

#### Start WGCNA #####

  phe <- pheno_cc[,-1]
  
  datExpr0 <- get(load(file= "/path/to/DE_downstream/WGCNA/input_VS_expression_filtered.Rdata"))
  datExpr0 <- t(datExpr0)
  pheno_cc <-pheno_cc[pheno_cc$rn %in% rownames(datExpr0), ]
  
  #check if genes have too many missings
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  if(!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  setwd("/path/to/DE_downstream/WGCNA/")
  
  sampleTree = hclust(dist(datExpr0), method = "average");
  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)

  datExpr = datExpr0
  
  #traitData
  allTraits = phe
  
  # Create pheno data matching the expression data
  Samples = rownames(datExpr);
  traitRows = match(Samples, allTraits$rn);
  datTraits = allTraits[traitRows, -1];
  rownames(datTraits) = allTraits$rn[traitRows] 
  datTraits$Cocaine_at_death <- as.numeric(as.factor(datTraits$Cocaine_at_death))
    
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  
  save(datExpr, datTraits, file = "VS_Expr_dataInput_processed_filtered.RData")
  
#### RNA-seq data ####

dir0 <- "/path/to/RNA"
  
  setwd(dir0)
  options(stringsAsFactors = FALSE)
  
  #read pheno
  pheno_r <- read.csv("/path/to/RNA/DE_results/pheno_coc_at_death.txt", sep=";")
  phe_orig <- read.csv2("/path/to/pheno/CUD_randomised_20220707_layout.csv")
  phe_orig <-phe_orig[,c(1,13,14)]
  
  pheno_all <- merge(pheno_r,phe_orig,by="Brain_ID")
  pheno_all$pH[is.na(pheno_all$pH)] <- 6.335
  
  #Remove PCA outliers
  pheno_all <- pheno_all[!(pheno_all$Brain_ID %in% c("ID1","ID2")),]
  
  # Add cibersort cell type estimates
  ciber <- read.csv("/path/to/DE_downstream/CIBERSORT/Markers_CIBERSORT_snRNA_VS.txt", sep=";")
  ciber <- ciber[,c(1:7)]
  colnames(ciber)[1] <- "rn"
  pheno_cc <- merge(pheno_all,ciber,by="rn")
  
  
#### Start WGCNA #####
  phe <- pheno_cc[,-2]
  datExpr0 <- get(load(file= "/path/to/DE_downstream/WGCNA/input_VS_expression_filtered.Rdata"))
  datExpr0 <- t(datExpr0)
  
  phe <- phe[phe$rn %in% rownames(datExpr0),]
  
  #check if genes have too many missings
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  if(!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  
  sampleTree = hclust(dist(datExpr0), method = "average");

  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  datExpr = datExpr0
  
  #traitData
  allTraits = phe[,-c(8,9)]
  
  # Create pheno data matching the expression data
  Samples = rownames(datExpr);
  traitRows = match(Samples, allTraits$rn);
  datTraits = allTraits[traitRows, -1];
  rownames(datTraits) = allTraits$rn[traitRows] 
  datTraits$Cocaine_at_death <- as.numeric(as.factor(datTraits$Cocaine_at_death))

  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  
  save(datExpr, datTraits, file = "VS_Expr_dataInput_processed_filtered.RData")
  
  
#### Protein dataset ####
dir0 <- "/path/to/protein/data_analysis_results_V1/"
  
  setwd(dir0)
  options(stringsAsFactors = FALSE)
  
  #read pheno
  pheno <- read.csv("/path/to/protein/data_analysis_results_V1/pheno_merged_with_proteomics.txt", sep=";")
  pheno_coc_at_death <- read.csv("/path/to/RNA/3_DE_results/pheno_coc_at_death.txt", sep=";")
  pheno_coc_at_death <- pheno_coc_at_death[,c("Brain_ID","Cocaine_at_death")]
  colnames(pheno_coc_at_death)[1] <-"SampleID" 
  pheno2 <- merge(pheno,pheno_coc_at_death,by="SampleID",all.x=T)
  pheno2$Cocaine_at_death[is.na(pheno2$Cocaine_at_death)] <- 1
  
  #Remove PCA outliers from pheno
  pheno2 <- pheno2[!(pheno2$ID %in% c("ID1","ID2")),]
  
  #### Start WGCNA ####
  
  datExpr0 <- t(read.csv("/path/to/protein/data_analysis_results_V1/normalized_data.txt", row.names=1, sep=";"))

  #Remove PCA outliers from the expression data
  datExpr0 <- datExpr0[!(rownames(datExpr0) %in% c("ID1","ID2")),]
  
  #check if genes have too many missings
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  if(!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  setwd("WGCNA/")
  
  sampleTree = hclust(dist(datExpr0), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  datExpr = datExpr0
  
  #traitData
  allTraits = pheno2
  
  # Create pheno data matching the expression data
  Samples = rownames(datExpr);
  traitRows = match(Samples, allTraits$ID);
  datTraits = allTraits[traitRows, -1];
  
  datTraits = datTraits[,c("CUD","Age","pH","PMI","batch","Cocaine_at_death")]
  rownames(datTraits) = allTraits$ID[traitRows] 
  datTraits$batch<- as.numeric(as.factor(datTraits$batch))
  datTraits$Cocaine_at_death<- as.numeric(as.factor(datTraits$Cocaine_at_death))

  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  
  save(datExpr, datTraits, file = "VS_Prot_dataInput_processed.RData")
  
  
