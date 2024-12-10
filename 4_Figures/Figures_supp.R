#### Code for generating the VS multi-omics supplementary figures ####
# The analysis code for specific R packages e.g. DESeq2, WGCNA, hdWGCNA and CellChat has been adapted based on the vignette documentation of the packages. Details can be found in the individual analysis scripts.
# Author: Eric Zillich 
# last modification: EZ 2024-12-04

###########################################################################
#### Supplementary Figure 1 - PCA of bulk-level datasets and CIBERSORT #### 
###########################################################################

#### S1A #### manually created 
#### S1B #### see script 3_PCA.R for code to generate the PCA plots
#### S1C ####

# Import CIBERSORT estimates
results<- read.csv("/path/to/RNA/DE_downstream/CIBERSORT/Markers_CIBERSORT_snRNA_VS.txt", sep=";")
rownames(results)<-results$Sample

# prepare results for plotting
pheno <- read.csv("/path/to/RNA/DE_results/pheno.txt", sep=";")
pheno <- pheno[,c("Brain_ID","rn","CUD")]
pheno <- pheno[order(pheno$CUD,decreasing=T),]
pheno$CUD[pheno$CUD == "1"] <- paste0("CUD_",c(1:19))
pheno$CUD[pheno$CUD == "0"] <- paste0("Ctrl_",c(1:21))

rownames(pheno) <-pheno$rn

res2 <- merge(pheno,results,by=0)
rownames(res2) <- res2$CUD

res3 <- res2[,c(6:11)]

perc_df <- data.frame(samples=rep(rownames(res3),each=6),celltypes=rep(colnames(res3),times=40),value=as.vector(t(res3)))
perc_df$samples <- factor(perc_df$samples, levels=pheno$CUD)

# Plot the fraction of celltypes as stacked barplot
per <- ggplot(perc_df, aes(fill=celltypes, y=value, x=samples)) + 
  geom_bar(position="fill", stat="identity")+theme(element_text(size=10))+scale_fill_manual(values=c("#F7AB64","#006960","#3A9BCC","#74BA59", "#70305A","#E8326D"))+ylab("celltype proportion")+xlab(NULL)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = per,filename = "/path/to/Figures/S1C.pdf", width = 10, height=4)

#### S1D ####
# Correlate the CIBERSROT cell type proportion estimation with the true snRNA celltype percentages in the N=16 samples where both information is available
# Import snRNA-seq data and extract cell type distribution
library(Seurat)
seurat<-readRDS("/path/to/snRNAseq/2_data_processed/2_harmony_object.rds")
seurat <- JoinLayers(seurat)
Idents(seurat) <- "orig.ident"

seurat$orig.ident[seurat$orig.ident =="GEX1"] <- "ID1"
seurat$orig.ident[seurat$orig.ident =="GEX2"] <- "ID2"
seurat$orig.ident[seurat$orig.ident =="GEX3"] <- "ID3"
seurat$orig.ident[seurat$orig.ident =="GEX4"] <- "ID4"
seurat$orig.ident[seurat$orig.ident =="GEX5"] <- "ID5"
seurat$orig.ident[seurat$orig.ident =="GEX6"] <- "ID6"
seurat$orig.ident[seurat$orig.ident =="GEX7"] <- "ID7"
seurat$orig.ident[seurat$orig.ident =="GEX8"] <- "ID8"
seurat$orig.ident[seurat$orig.ident =="GEX9"] <- "ID9"
seurat$orig.ident[seurat$orig.ident =="GEX10"] <- "ID10"
seurat$orig.ident[seurat$orig.ident =="GEX11"] <- "ID11"
seurat$orig.ident[seurat$orig.ident =="GEX12"] <- "ID12"
seurat$orig.ident[seurat$orig.ident =="GEX13"] <- "ID13"
seurat$orig.ident[seurat$orig.ident =="GEX14"] <- "ID14"
seurat$orig.ident[seurat$orig.ident =="GEX15"] <- "ID15"
seurat$orig.ident[seurat$orig.ident =="GEX16"] <- "ID16"

dist <- as.data.frame.matrix(table(seurat$orig.ident,seurat$celltypes))
dist <- dist[,colnames(dist) %in% c("Astrocyte","D1-MSN","D2-MSN","GABAergic-1","GABAergic-2","GABAergic-3","Microglia","Oligodendrocyte","OPC")]
dist$Inh_MSN <- rowSums(dist[,c("D1-MSN","D2-MSN")])
dist$Inh_GABA <- rowSums(dist[,c("GABAergic-1","GABAergic-2","GABAergic-3")])
dist <- dist[,c("Astrocyte","Inh_MSN","Inh_GABA","Microglia","Oligodendrocyte","OPC")]
dist$sum <- rowSums(dist)

for(i in rownames(dist)){
  for(j in colnames(dist)[c(1:6)]){
    dist[i,j] <- dist[i,j]/dist[i,"sum"]
  }
}

dist <- dist[,-7]

# subset CIBERSORT prediction in bulk RNA-seq data
res2_1 <-  res2[res2$Brain_ID %in% paste0(c("ID"),c(1:16)),c(2,5:10)]
rownames(res2_1) <- res2_1$Brain_ID
res2_1 <- res2_1[,-1]

res2_1<-res2_1[match(rownames(dist),rownames(res2_1)),match(colnames(dist),colnames(res2_1))]

# perform correlation analysis 
diag <- diag(cor(t(res2_1), t(dist), method = "pearson")) 
df_cor <- data.frame(Sample=names(diag),R=unname(diag))

# Plot the results
library(ggpubr)
p1 <- ggdotchart(df_cor, x = "Sample", y = "R", ggtheme = theme_pubr(base_size = 10), add = "segment", ylim = c(0, 1), ylab = "Corrleation CIBERSORT estimates - snRNAseq celltype proportion (Pearson R)")+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+xlab("snRNA seq samples (N=16)")+geom_hline(aes(yintercept = median(df_cor$R)),size=0.5,linetype="dashed")+ggtitle("Median R = 0.86")
ggsave(plot = p1,filename = "/path/to/Figures/S1D.pdf", width = 5, height=5)

##############################################################################
#### Supplementary Figure 2 - Transcriptome-proteome correlation analysis ####
##############################################################################

#### S2A ####
# Calculate the correlation between expr and prot
logTPM_int <- read.csv("/path/to/RNA/mRNA_prot_corr/logTPM.txt", row.names=1, sep="") # Import RNA data
datProt_int <- read.csv("/path/to/RNA/mRNA_prot_corr/logTMTint.txt", row.names=1, sep="") # Import protein data
diag <- diag(cor(logTPM_int, datProt_int, method = "pearson")) 

df_cor <- data.frame(Sample=names(diag),R=unname(diag))

p1 <- ggdotchart(df_cor, x = "Sample", y = "R", ggtheme = theme_pubr(base_size = 10), add = "segment", ylim = c(0, 1), ylab = "mRNA-protein correlation (Pearson R)") +theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+xlab("Samples (N=39)")
ggsave("/path/to/Figures/S2A.pdf",p1,width=6,height=8)

#### S2B ####
RNA_mean <- rowMeans(logTPM_int)
prot_mean <- rowMeans(datProt_int)
dat <- data.frame(RNA=unname(RNA_mean),Protein=unname(prot_mean),row.names = names(RNA_mean))

p2 <- ggscatter(dat, x = "RNA", y = "Protein", 
                size=0.1,color="gray10",
                cor.coef = TRUE, fullrange = T,cor.method = "pearson",
                xlab = "log2 TPM - RNA", ylab = "log2 TMT intensity - Protein",title="Correlation by gene")+theme(plot.title = element_text(face = "bold"))
ggsave("/path/to/Figures/S2B.pdf",p2,width=6,height=6)


#### S2C ####
# Top 40 are top 1% of expressed proteins - show the correlation for the highly expressed proteins only
dat_high_expr <- dat[order(dat$Protein,decreasing = T),]
dat_high_expr <- dat_high_expr[c(1:40),]
p3 <- ggscatter(dat_high_expr, x = "RNA", y = "Protein", 
                  size=0.1,color="gray10",label = rownames(dat_high_expr),font.label = c(8, "plain"),repel = T,
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "log2 TPM - RNA", ylab = "log2 TMT intensity - Protein",title="Highly expressed proteins")+theme(plot.title = element_text(face = "bold"))
ggsave("/path/to/Figures/S2C.pdf",p3,width=6,height=6)


#### S2D ####
# Separate by CUD and Ctrl to see if correlation is different in CUD and Ctrl, for phe_ordered see script 5b_correlation_analysis_transcriptome_proteome.R
CUD <- phe_ordered$CUD[match(colnames(logTPM_int),phe_ordered$Brain_ID)]
colnames(logTPM_int) <- CUD 
colnames(datProt_int) <- CUD 

dat_CUD <- logTPM_int[,colnames(logTPM_int) == 1]
RNA_mean_dat_CUD <- rowMeans(dat_CUD)
dat_Ctrl <- logTPM_int[,colnames(logTPM_int) == 0]
RNA_mean_dat_Ctrl <- rowMeans(dat_Ctrl)

dat_CUD2 <- datProt_int[,colnames(datProt_int) == 1]
Prot_mean_dat_CUD <- rowMeans(dat_CUD2)
dat_Ctrl2 <- datProt_int[,colnames(datProt_int) == 0]
Prot_mean_dat_Ctrl <- rowMeans(dat_Ctrl2)

dat_CUD_comb <- data.frame(RNA=unname(RNA_mean_dat_CUD),Protein=unname(Prot_mean_dat_CUD),row.names = names(RNA_mean_dat_CUD))
dat_Ctrl_comb <- data.frame(RNA=unname(RNA_mean_dat_Ctrl),Protein=unname(Prot_mean_dat_Ctrl),row.names = names(RNA_mean_dat_Ctrl))

dat_CUD_comb$status <- "CUD"
rownames(dat_CUD_comb) <- paste0("CUD_",rownames(dat_CUD_comb))
dat_Ctrl_comb$status <- "Ctrl"
rownames(dat_Ctrl_comb) <- paste0("Ctrl_",rownames(dat_Ctrl_comb))

dat_combined <- rbind(dat_CUD_comb,dat_Ctrl_comb)

# Plot correlation separated by CUD/Ctrl status
p4 <- ggscatter(dat_combined , x = "RNA", y = "Protein", 
                  size=0.1,color = "status",
                  palette = c("gray10","red"),
                  cor.coef = F, fullrange = T,cor.method = "pearson",
                  xlab = "log2 TPM - RNA", ylab = "log2 TMT intensity - Protein",title="Correlation by CUD status")+theme(plot.title = element_text(face = "bold"))+stat_cor(aes(color = status))+theme(legend.position = "none")
ggsave("/path/to/Figures/S2D.pdf",p4,width=6,height=6)


#### S2E ####
# Correlation coefficients for each gene where RNA and protein data is available
diag2 <- diag(cor(t(logTPM_int), t(datProt_int), method = "pearson"))

# Calculate median
med <- median(diag2)

# add quantiles 1% and 99% for gene annotation 
q <- quantile(diag2, probs = c(0.01,0.5,0.99))

# Extract genes at both ends of the distribution
low <- sort(diag2[diag2<unname(q[1])])
high <- sort(diag2[diag2>unname(q[3])])

p5 <- gghistogram(diag2, bins = 50, theme = theme_pubr(base_size = 10), fill = "gray",
                  xlab = "mRNA-protein correlation (Pearson R)",ylab = "Number of genes") +xlim(-1,1)+geom_vline(xintercept=med,col="red")+geom_vline(xintercept=unname(q[1]),col="red",linetype="dashed")+geom_vline(xintercept=unname(q[3]),col="red",linetype="dashed")
ggsave("/path/to/Figures/S2E.pdf",p5,width=6,height=8)

#### S2F ####
# GSEA on the gene-level correlation coefficients 
diag2 <- sort(diag2,decreasing = TRUE)

gse <- gseGO(
  diag2,
  OrgDb=org.Hs.eg.db,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")

gse<- pairwise_termsim(gse)
gse_pos <- gse@result[gse@result$NES > 0,]
gse_pos2 <- gse_pos[order(gse_pos$NES,decreasing=T),]
gse_pos2 <- gse_pos2[gse_pos2$p.adjust<0.05,c("Description","NES","p.adjust")]
gse_pos2$dir <- "positive"

gse_neg <- gse@result[gse@result$NES < 0,]
gse_neg2 <- gse_neg[order(gse_neg$NES,decreasing=F),]
gse_neg2 <- gse_neg2[gse_neg2$p.adjust<0.05,c("Description","NES","p.adjust")]
gse_neg2 <- gse_neg2[c(1:10),] # extract the 10 most significant pathways for plotting
gse_neg2$dir <- "negative"

gse <- rbind(gse_pos2,gse_neg2)
stars_gen <- function(x){
  stars <- c( "***", "**", "*","")
  vec <- c(0, 0.001, 0.01, 0.05,1.01)
  i <- findInterval(x, vec)
  stars[i]
}

gse$lab <- stars_gen(gse$p.adjust)


library(Seurat)
gsea_plot <- ggbarplot(gse, x = "Description", y = "NES",
                       fill = "dir",     
                       color = "black",           
                       sort.val = "desc",      
                       sort.by.groups = FALSE, 
                       lab.size=2)+ylim(-5,5)+NoLegend()+xlab(NULL)+scale_fill_manual(values=c("#268989","#E43F3F"))+ylab("GSEA set statistic")+ geom_text(aes(y=0.8*abs(NES)/NES,label = lab), vjust = 0.5,hjust=0,size=6)+coord_flip()+theme_minimal() + theme(text=element_text(size=20))

ggsave("/path/to/Figures/S2F.pdf",gsea_plot,width=12,height=8)

# Check for difference of correlation coefficients between CUD and Ctrl to calculate a delta R value indicating synchronized/desynchronized correlation patterns in CUD/Ctrl
# Correlation coefficients for each gene where RNA and protein data is available - split by CUD and Ctrl 
diag_CUD <- diag(cor(t(dat_CUD), t(dat_CUD2), method = "pearson"))
diag_Ctrl <- diag(cor(t(dat_Ctrl), t(dat_Ctrl2), method = "pearson"))

# Create results dataframe
dat_diff <- data.frame(CUD=unname(diag_CUD),Ctrl=unname(diag_Ctrl))
rownames(dat_diff)<- names(diag_CUD)
dat_diff$diff <- diag_CUD-diag_Ctrl # diff represents the delta R value
dat_diff <- dat_diff[order(dat_diff$diff,decreasing = T),]
write.table(dat_diff,"/path/to/results/dat_diff.txt",row.names = T, col.names = T, quote=F, sep=";")

###############################################################################
#### Supplementary Figure 3 - WGCNA full module-trait correlation heatmaps ####
###############################################################################

#### S3A #### generated during miRNA network construction from net.rds object (see script 6b_WGCNA_automatic_network_construction.R)
# Example code for generating the dendrogram plot:
mergedColors = labels2colors(net$colors)
png("VS_dendogram.png", width = 6, height = 4,units="in",res=600)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Correlation heatmap created using moduleTraitCor, datTraits, and textMatrix (see script 6c_WGCNA_relate_to_traits.R)

#### S3B ####  All miRNA modules
pdf("/path/to/Figures/S3B.pdf", width = 6.5, height= 6.5)
par(mar = c(6, 8.5, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor[,c(1,6,2:5)],
               xLabels = names(datTraits)[c(1,6,2:5)],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[,c(1,6,2:5)],
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#### S3C ####  All RNA modules
pdf("/path/to/Figures/S3C.pdf", width = 9, height = 19)
par(mar = c(6, 8.5, 5, 5))
labeledHeatmap(Matrix = moduleTraitCor[,c(1,6,2:5,7:12)],
               xLabels = names(datTraits)[c(1,6,2:5,7:12)],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[,c(1,6,2:5,7:12)],
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#### S3D ####  All protein modules
pdf("/path/to/Figures/S3D.pdf", width = 5, height = 8)
labeledHeatmap(Matrix = moduleTraitCor[,c(1,6,2:5)],
               xLabels = names(datTraits)[c(1,6,2:5)],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[,c(1,6,2:5)],
               setStdMargins = T,
               cex.text = 0.5,
               zlim = c(-1,1),
               cex.lab.y = 0.7,
               main = paste("Module-trait relationships"))
dev.off()

#############################################################
#### Supplementary Figure 4 - MOFA Supplementary Figures ####
#############################################################

# Import MOFA model (see script 7c_MOFA_downstream.Rmd)
VS <- readRDS("/path/to/protein/data_analysis_results_V1/MOFA/output/MOFA_VS_trained_model_20240303.rds")

# metadata
pheno_all_prot <- read.csv("/path/to/protein/data_analysis_results_V1/MOFA/pheno_all_prot.txt", sep=";")
pheno_MDD_AUD_coc_at_death <- read.csv("/path/to/RNA/DE_results/pheno.txt", sep=";")
pheno2 <- pheno_MDD_AUD_coc_at_death[,c("Brain_ID","Cocaine_at_death")]
pheno2[41,] <- c("ID3",1) 
pheno_all_prot <- merge(pheno_all_prot,pheno2,by="Brain_ID")

# Add cibersort cell type estimates
ciber <- read.csv("/path/to/RNA/DE_downstream/CIBERSORT/Markers_CIBERSORT_snRNA_VS.txt", sep=";")
ciber <- ciber[,c(1:7)]
colnames(ciber)[1] <- "rn_mRNA"

pheno_cc <- merge(pheno_all_prot,ciber,by="rn_mRNA",all.x=T)

colnames(pheno_cc)[2] <- "sample"
pheno_cc$sample <- as.character(pheno_cc$sample)
colnames(pheno_cc) <- gsub(".x$","",colnames(pheno_cc))
pheno_cc <-pheno_cc[,c(2,1,3:7,10:19)]
colnames(pheno_cc)[11] <- "Cocaine at death"
samples_metadata(VS) <- pheno_cc

#### S4A ####
library(corrplot)
pdf("/path/to/Figures/S4A.pdf", width = 5, height = 5)
plot_factor_cor(VS,col=rev(COL2("RdBu",200)))
dev.off()

#### S4B ####
var_exp <- plot_variance_explained(VS, plot_total = T)[[2]]
ggsave("/path/to/Figures/S4B.pdf",var_exp, width = 5, height = 4)

#### S4C ####
p_coc_at_death <- plot_factor(VS, 
                   factors = 10, 
                   color_by = "Cocaine_at_death",
                   add_violin = F,
                   dodge = TRUE
)+scale_fill_manual(values=c("#268989","#E43F3F"))+ylim(-2,3)
ggsave("/path/to/Figures/S4C.pdf",p_coc_at_death , width = 4, height = 3)

#### S4D ####
p_astro <- plot_factor(VS, 
                   factors = 10, 
                   color_by = "Astrocyte",
                   add_boxplot = TRUE,
                   dodge = TRUE
)+ylim(-2,3)
ggsave("/path/to/Figures/S4D.pdf",p_astro, width = 4, height = 3)

#### S4E ####
p_GABA <- plot_factor(VS, 
                   factors = 10, 
                   color_by = "Inh_GABA",
                   add_boxplot = TRUE,
                   dodge = TRUE
)+ylim(-2,3)
ggsave("/path/to/Figures/S4E.pdf",p_GABA, width = 4, height = 3)

#### S4F ####
p_oligo <- plot_factor(VS, 
                   factors = 10, 
                   color_by = "Oligodendrocyte",
                   add_boxplot = TRUE,
                   dodge = TRUE
)+ylim(-2,3)
ggsave("/path/to/Figures/S4F.pdf",p_oligo , width = 4, height = 3)

############################################################
#### Supplementary Figure 5 - Covariate plots snRNA-seq ####
############################################################

# Import Seurat object
seurat <- readRDS("/path/to/snRNAseq/2_data_processed/2_harmony_object.rds")
Idents(seurat) <- "celltypes"

library(viridis)
#### S5A ####
Idents(seurat)<- "CUD"
p_CUD <- DimPlot(seurat,label=F,reduction = "umap.harmony",label.size = 3)+xlab("UMAP1")+ylab("UMAP2")
ggsave("/path/to/Figures/S5A.pdf",plot = p_CUD,width = 6,height=5)

#### S5B ####
Idents(seurat)<- "orig.ident"
p_samp <- DimPlot(seurat,label=F,reduction = "umap.harmony",label.size = 3)+xlab("UMAP1")+ylab("UMAP2")
ggsave("/path/to/Figures/S5B.pdf",plot = p_samp,width = 6,height=5)

#### S5C ####
pheno_coc_at_death <- read.csv("/path/to/RNA/DE_results/pheno_coc_at_death.txt", sep=";")
seurat$coc_at_death <-""
seurat$coc_at_death[seurat$Brain_ID %in% pheno_coc_at_death$Brain_ID[pheno_coc_at_death$Cocaine_at_death==1]] <- "yes"
seurat$coc_at_death[seurat$Brain_ID %in% pheno_coc_at_death$Brain_ID[pheno_coc_at_death$Cocaine_at_death==0]] <- "no"

Idents(seurat)<- "coc_at_death"
p_coc_at_death2 <- DimPlot(seurat,label=F,reduction = "umap.harmony",label.size = 3)+xlab("UMAP1")+ylab("UMAP2")
ggsave("/path/to/Figures/S5C.pdf",plot = p_coc_at_death2 ,width = 6,height=5)

#### S5D ####
Idents(seurat)<- "Age"
seurat$Age <- as.character(seurat$Age)
seurat$Age <- factor(seurat$Age,levels=names(table(sort(as.numeric(seurat$Age)))))
p_age <- DimPlot(seurat,label=F,reduction = "umap.harmony",order=rev(levels(seurat$Age)))+xlab("UMAP1")+ylab("UMAP2")+scale_color_manual(values=rev(viridis(16)))
ggsave("/path/to/Figures/S5D.pdf",plot = p_age,width = 6,height=5)

#### S5E ####
Idents(seurat)<- "pH"
seurat$pH <- as.character(seurat$pH)
seurat$pH <- factor(seurat$pH,levels=names(table(sort(as.numeric(seurat$pH)))))
p_pH <- DimPlot(seurat,label=F,reduction = "umap.harmony",order=rev(levels(seurat$pH)))+xlab("UMAP1")+ylab("UMAP2")+scale_color_manual(values=rev(viridis(16)))
ggsave("/path/to/Figures/S5E.pdf",plot = p_pH,width = 6,height=5)

#### S5F ####
Idents(seurat)<- "PMI"
seurat$PMI <- as.character(seurat$PMI)
seurat$PMI <- factor(seurat$PMI,levels=names(table(sort(as.numeric(seurat$PMI)))))
p_pmi <- DimPlot(seurat,label=F,reduction = "umap.harmony",order=rev(levels(seurat$PMI)))+xlab("UMAP1")+ylab("UMAP2")+scale_color_manual(values=viridis(16))
ggsave("/path/to/Figures/S5F.pdf",plot = p_pmi,width = 6,height=5)

##############################################################################################################
#### Supplementary Figure 6 - Annotation of GABAergic clusters 1-3 and nuclei per condition per cell type ####
##############################################################################################################

#### S6A ####

Idents(seurat)<-"celltypes"

seurat_gaba <- subset(seurat, idents=c("GABAergic-1","GABAergic-2","GABAergic-3"))

# Expression marker gene heatmap
genes_gaba <- c("GAD1","GAD2","DRD1","ADARB2","CXCL14","CDH10","PTHLH","THSD4","NPY","SST","VIP","CCK","TAC3","PTPRK","TMEM163","GFRA2","CHAT","DRD2","HTR7","GRIK3","PTPRT")
Idents(seurat_gaba) <- "celltypes"
seurat_gaba@active.ident <- factor(seurat_gaba@active.ident,
                                   levels= rev(c("GABAergic-1",
                                                 "GABAergic-2",
                                                 "GABAergic-3")))

dp <- DotPlot(object = seurat_gaba, features = genes_gaba,assay ="RNA", dot.scale = 4.5) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=10)) + 
  theme(text=element_text(size = 8))+
  theme(aspect.ratio =0.8) + xlab("")+ylab("")
ggsave("/path/to/Figures/S6A.pdf",plot = dp,width = 7,height=5)

#### S6B ####
seurat$CUD <- factor(seurat$CUD,levels=c("CUD","Ctrl"))
library(gridExtra)
pdf("/path/to/Figures/S6B.pdf",width = 14,height=2)
grid.table(table(seurat$CUD,seurat$celltypes))
dev.off()

################################################################
#### Supplementary Figure 7 - snRNA-seq DE comparison plots ####
################################################################

# Import DE results - log2FC 0.5 cutoff and padj 0.001, only celltypes with more than 100 cells from each condition
OPC <- read.csv("/path/to/snRNAseq/3_results/DE/clusterOPC.csv")
OPC_sig <- OPC[OPC$p_val_adj<0.001 & abs(OPC$avg_log2FC)>0.5 ,]
OPC_sig <- OPC_sig[OPC_sig$pct.1 > 0.25 | OPC_sig$pct.2 > 0.25, ]

Oligodendrocyte <- read.csv("/path/to/snRNAseq/3_results/DE/clusterOligodendrocyte.csv")
Oligodendrocyte_sig <- Oligodendrocyte[Oligodendrocyte$p_val_adj<0.001 & abs(Oligodendrocyte$avg_log2FC)>0.5,]
Oligodendrocyte_sig <- Oligodendrocyte_sig[Oligodendrocyte_sig$pct.1 > 0.25 | Oligodendrocyte_sig$pct.2 > 0.25, ]

Astrocyte <- read.csv("/path/to/snRNAseq/3_results/DE/clusterAstrocyte.csv")
Astrocyte_sig <- Astrocyte[Astrocyte$p_val_adj<0.001 & abs(Astrocyte$avg_log2FC)>0.5,]
Astrocyte_sig <- Astrocyte_sig[Astrocyte_sig$pct.1 > 0.25 | Astrocyte_sig$pct.2 > 0.25, ]

Microglia <- read.csv("/path/to/snRNAseq/3_results/DE/clusterMicroglia.csv")
Microglia_sig <- Microglia[Microglia$p_val_adj<0.001 & abs(Microglia$avg_log2FC)>0.5,]
Microglia_sig <- Microglia_sig[Microglia_sig$pct.1 > 0.25 | Microglia_sig$pct.2 > 0.25, ]

GABAergic_1 <- read.csv("/path/to/snRNAseq/3_results/DE/clusterGABAergic-1.csv")
GABAergic_1_sig <- GABAergic_1[GABAergic_1$p_val_adj<0.001 & abs(GABAergic_1$avg_log2FC)>0.5,]
GABAergic_1_sig <- GABAergic_1_sig[GABAergic_1_sig$pct.1 > 0.25 | GABAergic_1_sig$pct.2 > 0.25, ]

D2MSN <- read.csv("/path/to/snRNAseq/3_results/DE/clusterD2-MSN.csv")
D2MSN_sig <- D2MSN[D2MSN$p_val_adj<0.001 & abs(D2MSN$avg_log2FC)>0.5,]
D2MSN_sig <- D2MSN_sig[D2MSN_sig$pct.1 > 0.25 | D2MSN_sig$pct.2 > 0.25, ]

D1MSN <- read.csv("/path/to/snRNAseq/3_results/DE/clusterD1-MSN.csv")
D1MSN_sig <- D1MSN[D1MSN$p_val_adj<0.001 & abs(D1MSN$avg_log2FC)>0.5,]
D1MSN_sig <- D1MSN_sig[D1MSN_sig$pct.1 > 0.25 | D1MSN_sig$pct.2 > 0.25, ]

#### S7A ####
#Upset plot for DE genes 
library(UpSetR)
cell_list <- list(Astrocyte = Astrocyte_sig$X,OPC = OPC_sig$X,Oligodendrocyte=Oligodendrocyte_sig$X,Microglia=Microglia_sig$X,D1MSN=D1MSN_sig$X,D2MSN=D2MSN_sig$X,GABAergic_1=GABAergic_1_sig$X)

pdf("/path/to/Figures/S7A.pdf",width=9,height=6.5)
upset(fromList(cell_list),nsets = 9,nintersects = 30,sets.x.label = "DE genes",mainbar.y.label = "DE gene overlap",order.by = "freq",text.scale = c(2,2,2,2,2,1.5))
dev.off()

#### S7B ####
# Run RRHO comparing snRNA-seq (cluster-ignorant) with bulk-level DE results
library(RRHO2)
markers <- FindMarkers(seurat,ident.1 = "CUD", ident.2 = "Ctrl", min.pct = 0.1,logfc.threshold = 0,assay = "RNA",features=genes) # perform DE analysis in snRNA-seq across all clusters (CUD vs. Ctrl comparison)

# Import bulk-level DE results
DE_VS <- read.csv("/path/to/RNA/DE_results/DE_results.txt")
DE_VS$RRHO_score <- -log10(DE_VS$pvalue)*sign(DE_VS$log2FoldChange)

# Prepare data for RRHO
markers$Gene <- rownames(markers)
markers$RRHO_score <- -log10(markers$p_val)*sign(markers$avg_log2FC)

genes_RRHO <- intersect(DE_VS$Gene,markers$Gene)

DE_VS <- DE_VS[DE_VS$Gene %in% 
                 genes_RRHO , c("Gene","RRHO_score")]
DE_snRNA_VS <- markers[markers$Gene %in% 
                         genes_RRHO , c("Gene","RRHO_score")]

# Run RRHO
RRHO_obj <-  RRHO2_initialize(DE_VS,DE_snRNA_VS , labels = c("bulk RNA-seq VS", "snRNA-seq VS"), log10.ind=TRUE,method = "hyper")

# Plot RRHO heatmap
pdf("/path/to/Figures/S7B.pdf",width=6,height=6)
RRHO2_heatmap(RRHO_obj)
dev.off()

#########################################################################
#### Supplementary Figure 8 - hdWGCNA in the human snRNA-seq dataset ####
#########################################################################

# Import hdWGCNA object 
seurat_obj <- readRDS("/path/to/snRNAseq/2_data_processed/3_hdWGCNA_object.rds")

#### S8A ####
# ModuleFeaturePlot for Astrocyte, Inh_MSN, and Inh_GABA 
p_Astrocyte <- ModuleFeaturePlot(seurat_obj, order='shuffle', raster=TRUE, raster_dpi=100, alpha=1, raster_scale=0.25, wgcna_name='Astrocyte', restrict_range=FALSE,reduction = "umap.harmony")

pdf(paste0(data_dir, "S8A_1.pdf"),height=20, width=20)
wrap_plots(p_Astrocyte, ncol=5)
dev.off()

p_Inh_GABA <- ModuleFeaturePlot(seurat_obj, order='shuffle', raster=TRUE, raster_dpi=100, alpha=1,raster_scale=0.25, wgcna_name='Inh_GABA', restrict_range=FALSE,reduction = "umap.harmony")

pdf(paste0(data_dir, "S8A_2.pdf"),height=20, width=20)
wrap_plots(p_Inh_GABA, ncol=5)
dev.off()

p_Inh_MSN <- ModuleFeaturePlot(seurat_obj, order='shuffle', raster=TRUE, raster_dpi=100, alpha=1,raster_scale=0.25, wgcna_name='Inh_MSN', restrict_range=FALSE,reduction = "umap.harmony")

pdf(paste0(data_dir, "S8A_3.pdf"),height=20, width=20)
wrap_plots(p_Inh_MSN, ncol=5)
dev.off()

#### S8B ####
# CUD DME analysis lollipop plot
dmel <- PlotDMEsLollipop2(
  seurat_obj, 
  CUD_DMEs[CUD_DMEs$group=="Astrocyte",], 
  wgcna_name="Astrocyte", 
  pvalue = "p_val_adj")+xlim(-1,2.5)

ggsave("/path/to/Figures/S8B.pdf",dmel, width=8, height=5)

#### S8C ####
dmel2 <- PlotDMEsLollipop2(
  seurat_obj, 
  CUD_DMEs[CUD_DMEs$group=="Inh_MSN",], 
  wgcna_name="Inh_MSN", 
  pvalue = "p_val_adj")+xlim(-1,2.8)
ggsave("/path/to/Figures/S8C.pdf",dmel2, width=8, height=5)

#### S8D ####
dmel3 <- PlotDMEsLollipop2(
  seurat_obj, 
  CUD_DMEs[CUD_DMEs$group=="Inh_GABA",], 
  wgcna_name="Inh_GABA", 
  pvalue = "p_val_adj")+xlim(-1,1)
ggsave("/path/to/Figures/S8D.pdf",dmel3, width=8, height=5)

#### S8E ####
# DE genes were loaded from Figure S7 (|log2FC|>0.5 and p.adj<0.001) and grouped based on positive (up) or negative (down) log2FC, now plot hdWGCNA module gene overlap with cluster-specific DE genes 

# Create DE gene list
l_DE <- list(Astrocyte_up_DE=Astrocyte_sig_up,Astrocyte_down_DE=Astrocyte_sig_down,D1MSN_up_DE=D1MSN_sig_up,D1MSN_down_DE=D1MSN_sig_down,D2MSN_up_DE=D2MSN_sig_up,D2MSN_down_DE=D2MSN_sig_down,Inh_GABA_up_DE=GABAergic_1_sig_up,Inh_GABA_down_DE=GABAergic_1_sig_down)

# Extract module genes for CUD-associated modules
# Astrocyte
modules <- GetModules(seurat_obj,wgcna_name = "Astrocyte")
g_ast_12 <- modules$gene_name[modules$module == "Astrocyte-M12"]
g_ast_14 <- modules$gene_name[modules$module == "Astrocyte-M14"]
g_ast_18 <- modules$gene_name[modules$module == "Astrocyte-M18"]

# Inh_MSN
modules <- GetModules(seurat_obj,wgcna_name = "Inh_MSN")
g_Inh_MSN_1 <- modules$gene_name[modules$module == "Inh_MSN-M1"]
g_Inh_MSN_2 <- modules$gene_name[modules$module == "Inh_MSN-M2"]
g_Inh_MSN_7 <- modules$gene_name[modules$module == "Inh_MSN-M7"]
g_Inh_MSN_10 <- modules$gene_name[modules$module == "Inh_MSN-M10"]

# Inh_GABA
modules <- GetModules(seurat_obj,wgcna_name = "Inh_GABA")
g_Inh_GABA_6 <- modules$gene_name[modules$module == "Inh_GABA-M6"]

# Create list of module genes
l_hdWGCNA <- list(Astocyte_M12=g_ast_12,Astrocyte_M14=g_ast_14,Astrocyte_M18=g_ast_18,Inh_MSN_M1=g_Inh_MSN_1,Inh_MSN_M2=g_Inh_MSN_2,Inh_MSN_M7=g_Inh_MSN_7,Inh_MSN_M10=g_Inh_MSN_10,Inh_GABA_M6=g_Inh_GABA_6)

# Run GeneOverlap analysis
library(GeneOverlap)
gom.obj <- newGOM(l_DE, l_hdWGCNA,36588)

pdf("/path/to/Figures/S8E.pdf",height=8, width=12)
drawHeatmap(gom.obj,grid.col="Reds", note.col="black",adj.p = T,what=c("Jaccard"))
dev.off()

#### S8F ####
# Hub gene network plot for Astrocyte, Inh_MSN, and Inh_GABA modules

groups <- c("Astrocyte","Inh_MSN","Inh_GABA")

for(cur_group in groups){
  ModuleNetworkPlot(
    seurat_obj,
    mods = "all",
    outdir = paste0(data_dir, cur_group, '_hubNetworks/'),
    wgcna_name = cur_group
  )
}

####################################################################################################
#### Supplementary Figure 9 - hdWGCNA GO/KEGG enrichment analysis results for the human dataset ####
####################################################################################################

# Run GO/KEGG enrichment analyses
dbs <-c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021', 'KEGG_2021_Human')

groups <- as.character(unique(seurat_obj$celltypes))
enrich_list <- list()

for(cur_group in groups){
  seurat_obj <- RunEnrichr(seurat_obj, dbs=dbs, wgcna_name=cur_group)
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name=cur_group) %>% subset(P.value < 0.05)
  enrichr_df$cell_type <- cur_group
  enrich_list[[cur_group]] <- enrichr_df
  
  EnrichrBarPlot(
    seurat_obj,
    outdir = paste0(data_dir,"enrichr_plots/",cur_group), 
    n_terms = 20, 
    plot_size = c(5,10),
    logscale=TRUE, wgcna_name=cur_group)
}


# Make GO/KEGG enrichment plots 

# Astrocyte

terms_ast<-GetEnrichrTable(seurat_obj,wgcna_name = "Astrocyte")
terms_ast_GO<-terms_ast[terms_ast$db %in% c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021"),]
terms_ast_KEGG<-terms_ast[terms_ast$db %in% c("KEGG_2021_Human"),]

#  Astrocyte-M12
m12 <- terms_ast_GO[terms_ast_GO$module=="Astrocyte-M12",]
m12_full <- m12[order(m12$Combined.Score,decreasing = T),][c(1:10),]
m12 <- m12[order(m12$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m12k <- terms_ast_KEGG[terms_ast_KEGG$module=="Astrocyte-M12",]
m12_k_full <- m12k[order(m12k$Combined.Score,decreasing = T),][c(1:10),]
m12k <- m12k[order(m12k$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m12 <- rbind(m12,m12k)
m12$Term[c(1:10)] <- gsub(".{13}$","",m12$Term[c(1:10)])
m12$db[c(1:10)]<-"GO"
m12$db[c(11:20)]<-"KEGG"

#  Astrocyte-M14
m14 <- terms_ast_GO[terms_ast_GO$module=="Astrocyte-M14",]
m14_full <- m14[order(m14$Combined.Score,decreasing = T),][c(1:10),]
m14 <- m14[order(m14$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m14k <- terms_ast_KEGG[terms_ast_KEGG$module=="Astrocyte-M14",]
m14_k_full <- m14k[order(m14k$Combined.Score,decreasing = T),][c(1:10),]
m14k <- m14k[order(m14k$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m14 <- rbind(m14,m14k)
m14$Term[c(1:10)] <- gsub(".{13}$","",m14$Term[c(1:10)])
m14$db[c(1:10)]<-"GO"
m14$db[c(11:20)]<-"KEGG"

# Inh_MSN
terms_msn<-GetEnrichrTable(seurat_obj,wgcna_name = "Inh_MSN")
terms_msn_GO<-terms_msn[terms_msn$db %in% c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021"),]
terms_msn_KEGG<-terms_msn[terms_msn$db %in% c("KEGG_2021_Human"),]

# Inh_MSN-M2
m2 <- terms_msn_GO[terms_msn_GO$module=="Inh_MSN-M2",]
m2_full <- m2[order(m2$Combined.Score,decreasing = T),][c(1:10),]
m2 <- m2[order(m2$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m2k <- terms_msn_KEGG[terms_msn_KEGG$module=="Inh_MSN-M2",]
m2_k_full <- m2k[order(m2k$Combined.Score,decreasing = T),][c(1:10),]
m2k <- m2k[order(m2k$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m2 <- rbind(m2,m2k)
m2$Term[c(1:10)] <- gsub(".{13}$","",m2$Term[c(1:10)])
m2$db[c(1:10)]<-"GO"
m2$db[c(11:20)]<-"KEGG"

# Inh_MSN-M7
m7 <- terms_msn_GO[terms_msn_GO$module=="Inh_MSN-M7",]
m7_full <- m7[order(m7$Combined.Score,decreasing = T),][c(1:10),]
m7 <- m7[order(m7$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m7k <- terms_msn_KEGG[terms_msn_KEGG$module=="Inh_MSN-M7",]
m7_k_full <- m7k[order(m7k$Combined.Score,decreasing = T),][c(1:10),]
m7k <- m7k[order(m7k$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m7 <- rbind(m7,m7k)
m7$Term[c(1:10)] <- gsub(".{13}$","",m7$Term[c(1:10)])
m7$db[c(1:10)]<-"GO"
m7$db[c(11:20)]<-"KEGG"

# Generate dot plots for the enrichment results

library(cowplot)

#### S9A ####
#Astrocyte-M12
m12$Combined.Score <- log(m12$Combined.Score)
position1 <- rev(m12$Term)

p1 <- ggplot(data = m12[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="red") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(2,7)
p2 <- ggplot(data = m12[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="red") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(2,7)

pdf("/path/to/Figures/S9A.pdf",width=10,height=9)
plot_grid(p1,p2,ncol=1,align="v")
dev.off()

#### S9B ####
#Astrocyte-M14
m14$Combined.Score <- log(m14$Combined.Score)
position1 <- rev(m14$Term)

p3 <- ggplot(data = m14[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="green") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(2,8)
p4 <- ggplot(data = m14[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="green") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(2,8)

pdf("/path/to/Figures/S9B.pdf",width=10,height=9)
plot_grid(p3,p4,ncol=1,align="v")
dev.off()

#### S9C ####
#Inh_MSN-M2
m2$Combined.Score <- log(m2$Combined.Score)
position1 <- rev(m2$Term)

p5 <- ggplot(data = m2[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="turquoise") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(1,6)
p6 <- ggplot(data = m2[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="turquoise") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(1,6)

pdf("/path/to/Figures/S9C.pdf",width=10,height=9)
plot_grid(p5,p6,ncol=1,align="v")
dev.off()

#### S9D ####
#Inh_MSN-M7
m7$Combined.Score <- log(m7$Combined.Score)
position1 <- rev(m7$Term)

p7 <- ggplot(data = m7[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="#FEE12B") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(2,7)
p8 <- ggplot(data = m7[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="#FEE12B") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(2,7)

pdf("/path/to/Figures/S9D.pdf",width=10,height=9)
plot_grid(p7,p8,ncol=1,align="v")
dev.off()
