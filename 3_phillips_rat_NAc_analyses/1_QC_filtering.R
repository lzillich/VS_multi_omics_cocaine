#### QC and filtering of Phillips et al. snRNA-seq data (obtained from NCBI GEO,PMID: 36965548) ####
# Author: Eric Zillich
# Last modification: EZ 2024-09-13 

library(Seurat)
library(dplyr)
set.seed(42)

# QC 
for (i in c("Cocaine","Ctrl")){

counts <- Read10X(paste0("/path/to/Phillips_VS_data/1_processed_GEO/",i,"/"))

# Create Seurat object containing RNA and ATAC data
sample <- CreateSeuratObject(counts = counts)
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^Mt-")

# Removal of low quality cells
sample <- subset(sample,
                 subset = nFeature_RNA > 900 &
                   nFeature_RNA < 8500 &
                   percent.mt < 10)


saveRDS(sample,file = paste0("/path/to/Phillips_VS_data/2_data_processed/postQC/",i,".rds"))
}


