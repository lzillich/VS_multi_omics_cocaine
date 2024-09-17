#### PCA for miRNA-seq, RNA-seq, and proteomics datasets #### 
# Author: Eric Zillich
# Last change 2023-07-12 EP

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

#### miRNA-seq ####
#### 1. Import phenotype data ####
seq_info <- read.delim("/path/to/Sample Allocation P2023-034-miRNA.txt", header=FALSE)
colnames(seq_info) <- c("rn","Brain_ID")
pheno <- read.csv2("/path/to/pheno/CUD_randomised_20220707_layout.csv")
RIN <- read_excel("/path/to/DE_results/RNA_info_RIN.xlsx")
RIN <- RIN[c(1:40), c("ID","RIN")]
colnames(RIN) <- c("Brain_ID","RIN")
RIN$RIN <- as.numeric(RIN$RIN)

pheno2 <- merge(seq_info,pheno,by="Brain_ID")
pheno2 <- merge(pheno2,RIN,by="Brain_ID")

# Remove irrelevant pheno data 
pheno2 <- pheno2[,c("Brain_ID","rn","CUD","Age","PMI","pH","RIN")]
write.table(pheno2,"/path/to/DE_results/pheno.txt",sep=";",quote = F,row.names = T,col.names = T)

#mean impute missing pH data 
mean(pheno2$pH,na.rm=T)
pheno2$pH[is.na(pheno2$pH)] <- 6.335
pheno2$CUD[pheno2$CUD == "1"] <- "CUD"
pheno2$CUD[pheno2$CUD == "0"] <- "Ctrl"

DT <- data.table
DF <- data.frame

#### 2. Import counts ####
# CUD
miRNA_CUD <- get(load("/path/to/miRNA/2_quantified/miRNA_CUD_counts.Rdata"))
colnames(miRNA_CUD$counts) <- gsub("_1_trimmedAligned.sortedByCoord.out.bam","",colnames(miRNA_CUD$counts))
miRNA_CUD$targets <- gsub("_1_trimmedAligned.sortedByCoord.out.bam","",miRNA_CUD$targets)
colnames(miRNA_CUD$stat) <- gsub("_1_trimmedAligned.sortedByCoord.out.bam","",colnames(miRNA_CUD$stat))

miRNA_CUD_counts <- DF(miRNA_CUD$counts)

pheno_with_CUD_Ctrl <- read.csv("/path/to/RNA/DE_results/pheno_with_CUD_Ctrl.txt", sep=";")
colnames(miRNA_CUD_counts) <- pheno_with_CUD_Ctrl$CUD[match(colnames(miRNA_CUD_counts),pheno_with_CUD_Ctrl$rn)]

#### 3. START ANALYSIS ####
#pheno
pheno_miRNA <- merge(pheno2,pheno_with_CUD_Ctrl,by="Brain_ID")
pheno_miRNA <- pheno_miRNA[order(pheno_miRNA$CUD.y,decreasing = T),]
miRNA_CUD_counts <- miRNA_CUD_counts[,order(colnames(miRNA_CUD_counts),decreasing = T)]

pheno_miRNA$CUD.y  == colnames(miRNA_CUD_counts)
# TRUE, continue 

#coldata
coldata <- as.data.frame(pheno_miRNA)
coldata$CUD <- as.factor(coldata$CUD.x)
coldata$Age <- scale(coldata$Age)
coldata$PMI <- scale(coldata$PMI)
coldata$RIN <- scale(coldata$RIN)
coldata$pH <- scale(coldata$pH)

rownames(coldata) <- coldata$CUD.y

#sort
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

## Define contrasts, extract results table and shrink log2 fold changes
contrast <-  c("CUD", "1", "0")

res <- results(dds, contrast=contrast, alpha = 0.05)

#### 4. Run PCA to get a broad overview of the dataset ####
vst <- vst(dds,nsub=828) #First do variance stabilizing transformation 
rv <- rowVars(assay(vst))

# select top variable genes
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform PCA
pca <- prcomp(t(assay(vst)[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
intgroup = "CUD.x"
intgroup.df <- as.data.frame(colData(vst)[, intgroup, drop=FALSE])

group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(vst)[[intgroup]]
}

# assemble the data frame for plotting
d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(vst))

pdf("/path/to/miRNA/DE_results/PCA_miRNA_Suppl.pdf", width = 5, height = 5 )
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=1) + 
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylim(-60,60)+xlim(-60,60)+
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed(ratio=1)+ geom_text_repel(size = 3, aes(label = name), vjust = 2.2,max.overlaps = 3) +theme_minimal()+scale_color_manual(values=c("black","black"))
dev.off()


#### RNA-seq ####
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

pheno_with_CUD_Ctrl <- read.csv("/path/to/RNA/DE_results/pheno_with_CUD_Ctrl.txt", sep=";")
colnames(mRNA_CUD_counts) <- pheno_with_CUD_Ctrl$CUD[match(colnames(mRNA_CUD_counts),pheno_with_CUD_Ctrl$rn)]

#### 3. START ANALYSIS ####

#pheno
pheno_mRNA <- merge(pheno2,pheno_with_CUD_Ctrl,by="Brain_ID")
colnames(pheno_mRNA)[colnames(pheno_mRNA) == "CUD.y"] <- "CUD" 
pheno_mRNA <- pheno_mRNA[order(pheno_mRNA$CUD,decreasing = T),]
mRNA_CUD_counts <- mRNA_CUD_counts[,order(colnames(mRNA_CUD_counts),decreasing = T)]

pheno_mRNA$CUD  == colnames(mRNA_CUD_counts)
# TRUE, continue 

#coldata
coldata <- as.data.frame(pheno_mRNA)
coldata$CUD2 <- as.factor(coldata$CUD.x)
coldata$Age <- scale(coldata$Age)
coldata$PMI <- scale(coldata$PMI)
coldata$RIN <- scale(coldata$RIN)
coldata$pH <- scale(coldata$pH)

rownames(coldata) <- coldata$CUD
  
#sort
coldata <- coldata[colnames(mRNA_CUD_counts),]
identical(colnames(mRNA_CUD_counts),rownames(coldata)) # must be TRUE

## Specify model
dds <- DESeqDataSetFromMatrix(countData = mRNA_CUD_counts,
                              colData = coldata,
                              design = ~ CUD2 + Age + RIN + pH + PMI)
# Prefiltering
keep <- rowSums(counts(dds) >= 4) >=8
dds <- dds[keep,]

# DE analysis
dds <- DESeq(dds)

## Define contrasts, extract results table and shrink log2 fold changes
contrast <-  c("CUD2", "CUD", "Ctrl")
res <- results(dds, contrast=contrast, alpha = 0.05)

#### 4. Run PCA to get a broad overview of the dataset ####
vst <- vst(dds) # First do variance stabilizing transformation 
rv <- rowVars(assay(vst))

# select top variable genes
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform PCA 
pca <- prcomp(t(assay(vst)[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

intgroup = "CUD.x"
intgroup.df <- as.data.frame(colData(vst)[, intgroup, drop=FALSE])

group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(vst)[[intgroup]]
}

# assemble the data frame for plotting
d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(vst))

pdf("/path/to/RNA/DE_results/PCA_RNA_Suppl.pdf", width = 5, height = 5 )
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=1) + 
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylim(-60,60)+xlim(-60,60)+
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed(ratio=1)+ geom_text_repel(size = 3, aes(label = name), vjust = 2.2,max.overlaps = 3) +theme_minimal()+scale_color_manual(values=c("black","black"))
dev.off()


#### Proteomics ####
# For objects raw_dataE, norm_dataE, and ctrl.ratio_dataE, see 4_Proteomics.Rmd

sets <- list("raw data" = raw_dataE, "batchcl data" = batchcl_dataE, "norm data" = norm_dataE, "ctrl.ratio data" = ctrl.ratio_dataE)
PCA_data <- NULL
for (i in seq_along(sets)) {
  set <- sets[[i]]
  set.name <- names(sets)[i]
  PCA_m <- t(na.omit(exprs(set)))
  if (length(PCA_m) > 50) {
    PCA <- prcomp(PCA_m, scale = FALSE)
    perc_var <-
      round(100 * PCA$sdev ^ 2 /
              sum(PCA$sdev ^ 2), 1)
    PCA_data_i <-
      data.frame(PC1 = PCA$x[, 1],
                 PC2 = PCA$x[, 2],
                 PC3 = PCA$x[, 3],
                 PC1.var = perc_var[1],
                 PC2.var = perc_var[2],
                 condition = pData(set)$condition,
                 rep = pData(set)$batch,
                 measurement = set.name)
    PCA_data <- bind_rows(PCA_data, PCA_data_i)
    rm(PCA_data_i)
  }
}
rm(i)
PCA_data$measurement <- factor(PCA_data$measurement, ordered = TRUE, levels = c("raw data", "batchcl data", "norm data", "ctrl.ratio data"))

ggplot(data = PCA_data, aes(PC1, PC2, colour = PC3)) +
  geom_point(aes(colour = condition, shape = rep), size = 4) +
  guides(size = "none") +
  customPlot +
  ggtitle("PCA analysis") +
  xlab("PC1") +
  ylab("PC2") +
  facet_wrap(. ~ measurement + paste("PC1:", PC1.var, "% var - PC2:", PC2.var, "% var"), scale = "free", ncol = 5)
ggsave(file.path(dir_save, paste0("PCA_analysis_", script.version, ".pdf")),
       width = 17, height = 5)
write.csv(PCA_data, file = file.path(dir_save, paste0("PCA_analysis_data_", script.version, ".csv")), row.names = FALSE)

# PCA plot

PCA_data2 <- PCA_data[PCA_data$measurement =="norm data",]
rnames <- substr(rownames(PCA_data2),10,12)
rownames(PCA_data2) <- rnames

# Import additional pheno information
pheno_with_CUD_Ctrl <- read.csv("/path/to/RNA/DE_results/pheno_with_CUD_Ctrl.txt", sep=";")
rownames(pheno_with_CUD_Ctrl) <- pheno_with_CUD_Ctrl$Brain_ID

PCA_data2 <- merge(pheno_with_CUD_Ctrl,PCA_data2,by=0)

# Plot PCA
library(ggrepel)
pdf("/path/to/protein/data_analysis_results_V1/PCA_Proteomics_Supp.pdf", width = 5, height = 5 )
ggplot(data=PCA_data2, aes_string(x="PC1", y="PC2", color="condition")) + geom_point(size=1) + 
  xlab(paste0("PC1: ",round(PCA_data2$PC1.var[1]),"% variance")) +ylim(-60,60)+xlim(-60,60)+
  ylab(paste0("PC2: ",round(PCA_data2$PC2.var[1]),"% variance")) +
  coord_fixed(ratio=1)+ geom_text_repel(size = 3, aes(label = CUD), vjust = 2.2,max.overlaps = 3) +theme_minimal()+scale_color_manual(values=c("black","black"))
dev.off()