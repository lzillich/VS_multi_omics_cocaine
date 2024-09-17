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

DF <- data.frame

#### miRNA data ####
dir0 <- "/path/to/miRNA"

setwd(dir0)
options(stringsAsFactors = FALSE)

#read pheno
pheno <- read.csv("/path/to/miRNA/DE_results/pheno.txt", sep=";")
pheno_c<- read.csv("/path/to/RNA/DE_results/pheno_coc_at_death.txt", sep=";")
pheno_c<-pheno_c[,c(1,8)]

pheno$pH[is.na(pheno$pH)] <- 6.335
pheno_cc <- merge(pheno,pheno_c,by="Brain_ID")

# remove PCA outliers 
pheno_cc <- pheno_cc[!(pheno_cc$rn %in% c("DE80NGSUKBR136473","DE42NGSUKBR136478")),]

# Allow multi-threading within WGCNA.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads(nThreads = 10)


###### Run WGNCA on miRNA data ########

        setwd("DE_downstream/WGCNA/")
        lnames = load(file ="VS_Expr_dataInput_processed_filtered.RData")    
        # Choose a set of soft-thresholding powers
        powers = c(c(1:10), seq(from = 12, to=30, by=2))
        
        # Call the network topology analysis function
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        
        sft$powerEstimate
        
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
             xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
        text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
             labels = powers, col = "red")
        # this line corresponds to using an R^2 cut-off of h
        abline(h = 0.9, col = "red")
        # Mean connectivity as a function of different powers
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
             ylab = "Mean Connectivity", main = paste("Mean connectivity"))
        text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
        # 6 seems to be the best here
        
        net = blockwiseModules(datExpr, power = 6,
                               TOMType = "signed", minModuleSize = 10, maxBlockSize = 36000,
                               reassignThreshold = 0, mergeCutHeight = 0.15,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "/path/to/miRNA/VS_miRNA_TOM",
                               verbose = 3)
        table(net$colors)
        
        setwd("/path/to/miRNA/")
        
        # Convert labels to colors for plotting
        mergedColors = labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        png("VS_dendogram_Expr.png", width = 6, height = 4,units="in",res=600)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        dev.off()
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        save(MEs, moduleLabels, moduleColors, geneTree,
             file = "VS_networkConstruction_Expr.RData")
        
#### RNA data ####
dir0 <- "/path/to/RNA"
        
        setwd(dir0)
        options(stringsAsFactors = FALSE)
        enableWGCNAThreads(nThreads = 20)
        
        ###### Run WGCNA on RNA data ########
        
        setwd("DE_downstream/WGCNA/")
        lnames = load(file ="VS_Expr_dataInput_processed_filtered.RData")    
        # Choose a set of soft-thresholding powers
        powers = c(1:40)
        
        # Call the network topology analysis function
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        
        sft$powerEstimate # estimate is 4 but connectivity should be lower (around 100, so check manually)
        
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
             xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
        text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
             labels = powers, col = "red")
        # this line corresponds to using an R^2 cut-off of h
        abline(h = 0.9, col = "red")
        # Mean connectivity as a function of different powers
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
             ylab = "Mean Connectivity", main = paste("Mean connectivity"))
        text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
        
        # 7 seems to be the best here
        
        net = blockwiseModules(datExpr, power = 7,
                               TOMType = "signed", minModuleSize = 20, maxBlockSize = 36000,
                               reassignThreshold = 0, mergeCutHeight = 0.15,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "/path/to/RNA/VS_mRNA_TOM",
                               verbose = 3)
        saveRDS(net,"/path/to/RNA/net.rds")
        table(net$colors)
        
        setwd("/path/to/RNA/")
        
        # Convert labels to colors for plotting
        mergedColors = labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        png("VS_dendogram_Expr.png", width = 6, height = 4,units="in",res=600)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        dev.off()
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        save(MEs, moduleLabels, moduleColors, geneTree,
             file = "VS_networkConstruction_Expr.RData")
        
        
#### Protein data ####
        dir0 <- "/path/to/protein/data_analysis_results_V1/"
        
        setwd(dir0)
        options(stringsAsFactors = FALSE)
        
        #read pheno
        pheno <- read.csv("/path/to/protein/data_analysis_results_V1/pheno_merged_with_proteomics.txt", sep=";")
        pheno_coc_at_death <- read.csv("/path/to/RNA/DE_results/pheno_coc_at_death.txt", sep=";")
        pheno_coc_at_death <- pheno_coc_at_death[,c("Brain_ID","Cocaine_at_death")]
        colnames(pheno_coc_at_death)[1] <-"SampleID" 
        pheno2 <- merge(pheno,pheno_coc_at_death,by="SampleID",all.x=T)
        #214 also had cocaine at death
        pheno2$Cocaine_at_death[is.na(pheno2$Cocaine_at_death)] <- 1
        enableWGCNAThreads(nThreads = 20)
        
        ###### Run WGCNA on protein data ########
        
        setwd("WGCNA/")
        lnames = load(file ="VS_Prot_dataInput_processed.RData")    
        # Choose a set of soft-thresholding powers
        powers = c(1:30)
        
        # Call the network topology analysis function
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        
        sft$powerEstimate #7
        
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
             xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
        text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
             labels = powers, col = "red")
        # this line corresponds to using an R^2 cut-off of h
        abline(h = 0.9, col = "red")
        # Mean connectivity as a function of different powers
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
             ylab = "Mean Connectivity", main = paste("Mean connectivity"))
        text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
        
        
        net = blockwiseModules(datExpr, power = 9,
                               TOMType = "signed", minModuleSize = 20, maxBlockSize = 36000,
                               reassignThreshold = 0, mergeCutHeight = 0.15,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = T,
                               saveTOMFileBase = "/path/to/protein/VS_Protein_TOM",
                               verbose = 3)
        table(net$colors)
        
        setwd("/path/to/protein/")
        
        # Convert labels to colors for plotting
        mergedColors = labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        png("VS_dendogram_Prot.png", width = 6, height = 4,units="in",res=600)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        dev.off()
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        save(MEs, moduleLabels, moduleColors, geneTree,
             file = "VS_networkConstruction_Prot.RData")
        
        
        
