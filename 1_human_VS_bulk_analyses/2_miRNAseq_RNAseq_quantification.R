### Quantification CUD miRNA-seq/RNA-seq data ###
# Author: Eric Zillich
# Last modification: 2024-09-10 EZ

library(readr)
library(data.table)
library(dplyr)
library(varhandle)
library(Rsubread)

#### miRNA-seq ####
# Read in miRNA samples
CUD_miRNA_files <- list.files("/path/to/miRNA/STAR_aligned/",pattern="*.Coord.out.bam$")# Quantification miRNA Seq

count <- featureCounts(c(paste0("/path/to/miRNA/STAR_aligned/",CUD_miRNA_files)), annot.ext="/path/to/reference/hsa.gff3",isGTFAnnotationFile = T,GTF.featureType = "miRNA",GTF.attrType = "Name",isPairedEnd = F, nthreads = 10,countMultiMappingReads = T,fraction=T)

#round the counts that were created from fractional mapping of multipmapping reads 
count$counts <- apply(count$counts,2,round)
save(count, file = "/path/to/miRNA/quantified/miRNA_CUD_counts.Rdata")

#### RNA-seq ####
# Read in RNA samples
CUD_mRNA_files <- list.files("/path/to/RNA/1_STAR_aligned/",pattern="*.Coord.out.bam$")# Quantification mRNA Seq

count <- featureCounts(c(paste0("/path/to/RNA/1_STAR_aligned/",CUD_mRNA_files)), annot.ext="/path/to/reference/GRCh38.subset.gtf",isGTFAnnotationFile = T,GTF.featureType = "exon",GTF.attrType = "gene_name",isPairedEnd = T, nthreads = 10,countMultiMappingReads = F,reportReads = "CORE")
save(count, file = "/path/to/RNA/2_quantified/mRNA_CUD_counts.Rdata")



