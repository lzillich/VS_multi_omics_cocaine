#### QC and filtering of snRNA-seq data - SoupX code adapted from its vignette: https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html ####
# Author: Eric Zillich
# Last modification: EZ 2024-09-13 

library(Seurat)
library(dplyr)
library(SoupX)
set.seed(42)

# QC and SoupX removal of ambient RNA

for (i in c(paste0("GEX",1:16))){

counts <- Read10X(paste0("/path/to/snRNAseq/1_cellranger_out/",i,"/outs/filtered_feature_bc_matrix/"))

# Create Seurat object containing RNA and ATAC data
sample <- CreateSeuratObject(counts = counts)
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")

# Removal of low quality cells
sample <- subset(sample,
                   subset = nFeature_RNA > 900 &
                     nFeature_RNA < 8500 &
                     percent.mt < 10)

# SoupX correction for ambient RNA
raw <- Read10X(paste0("/path/to/snRNAseq/1_cellranger_out/",i,"/outs/raw_feature_bc_matrix/"))

soup.channel  <- SoupChannel(raw, sample[["RNA"]]$counts)
srat    <- sample
srat    <- NormalizeData(srat,verbose=F)
srat    <- ScaleData(srat,verbose = F)
srat    <- FindVariableFeatures(srat,verbose=F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat,dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel <- setClusters(soup.channel,setNames(meta$seurat_clusters,rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
soup.channel  <- autoEstCont(soup.channel)
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

sample[["origcounts"]] <- CreateAssayObject(counts=sample[["RNA"]]$counts)
sample[["RNA"]]$counts <- adj.matrix

saveRDS(sample,file = paste0("/path/to/snRNAseq/2_data_processed/postQC/",i,".rds"))
}


