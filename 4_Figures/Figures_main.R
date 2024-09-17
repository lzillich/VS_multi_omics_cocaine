#### Code for generating the VS multi-omics main manuscript figures ####
# The analysis code for specific R packages e.g. DESeq2, WGCNA, hdWGCNA and CellChat has been adapted based on the vignette documentation of the packages. Details can be found in the individual analysis scripts.
# Author: Eric Zillich 
# last modification: EZ 2024-09-13

#############################################################
#### Figure 1 - DE and GO analyses miRNA, mRNA, proteins ####
#############################################################

#### 1A+B #### generated using BioRender (https://www.biorender.com/)

#### 1C ####
library(ggrepel)
library(ggplot2)

DE_results <- read.csv("/path/to/miRNA/DE_results/DE_results.txt", row.names=1)
DE_results$miR <- rownames(DE_results)
DE_results$DE <- 
  with(DE_results,ifelse(pvalue < 0.05 & log2FoldChange > 0.07, "DE up p<0.05",
                         ifelse(pvalue < 0.05 & log2FoldChange < -0.07, "DE down p<0.05", "n.s.")))
DE_results$DE <- 
  factor(DE_results$DE, 
         ordered = TRUE, 
         levels = c("DE up p<0.05","DE down p<0.05","n.s."))

DE_results <- DE_results[order(DE_results$pvalue,decreasing = F),]
  

customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","gray80")), 
  scale_colour_manual(values=c("#E43F3F","#268989","gray80"))
)


p1 <- ggplot(data = DE_results, aes(log2FoldChange, -log10(pvalue), colour = DE)) + xlim(-6,6)+ ylim(0,6.7)+
  geom_point(size=0.2) + geom_hline(aes(yintercept = -log10(0.05)),size=0.2) + geom_hline(aes(yintercept = -log10(3.405995e-05)),size=0.2,linetype="dashed") + geom_vline(aes(xintercept=0),size=0.2)+
  geom_text_repel(aes(label = c(DE_results$miR[DE_results$log2FoldChange>0][c(1:2)],DE_results$miR[DE_results$log2FoldChange<0][c(1:2)]),segment.size=0.2), 
                  data = rbind(DE_results[DE_results$log2FoldChange>0,][c(1:2),],DE_results[DE_results$log2FoldChange<0,][c(1:2),]), 
                  vjust = 0, nudge_y = 0.1, size = 3) +
  xlab("log2FoldChange") +ylab("-log10(p)")+
  customPlot+theme(legend.position = "none")+ggtitle("microRNA")+theme(plot.title = element_text(hjust = 0.5))


#### 1D ####

DE_res <- read.csv("/path/to/RNA/DE_results/DE_results.txt")
DE_res$DE <- 
  with(DE_res,ifelse(pvalue < 0.05 & log2FoldChange > 0.07 , "DE up p<0.05",
                         ifelse(pvalue < 0.05 & log2FoldChange < -0.07 , "DE down p<0.05", "n.s.")))
DE_res$DE <- 
  factor(DE_res$DE, 
         ordered = TRUE, 
         levels = c("DE up p<0.05","DE down p<0.05","n.s."))

DE_res <- DE_res[order(DE_res$pvalue,decreasing = F),]
library(ggplot2)

customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","gray80")), 
  scale_colour_manual(values=c("#E43F3F","#268989","gray80"))
)


p2 <- ggplot(data = DE_res, aes(log2FoldChange, -log10(pvalue), colour = DE)) + xlim(-3,3)+ ylim(0,6.7)+
  geom_point(size=0.2) + geom_hline(aes(yintercept = -log10(0.05)),size=0.2) + geom_hline(aes(yintercept = -log10( 1.017346e-05)),size=0.2,linetype="dashed") + geom_vline(aes(xintercept=0),size=0.2)+
  geom_text_repel(aes(label = c(DE_res$Gene[DE_res$log2FoldChange>0][c(1:5)],DE_res$Gene[DE_res$log2FoldChange<0][c(1:5)]),segment.size=0.2), 
                  data = rbind(DE_res[DE_res$log2FoldChange>0,][c(1:5),],DE_res[DE_res$log2FoldChange<0,][c(1:5),]), 
                  vjust = 0, nudge_y = 0.1, size = 3) +
  xlab("log2FoldChange") +ylab("")+
  customPlot+theme(legend.position = "none")+ggtitle("mRNA")+theme(plot.title = element_text(hjust = 0.5))


#### 1E ####
limma_results <- read.csv("/path/to/protein/data_analysis_results_V1/Limma_results_V1.csv")

limma_results$DE <- 
  with(limma_results, ifelse(pvalue.limma < 0.05 & logFC > 0.07, "DE up p<0.05",
                             ifelse(pvalue.limma < 0.05 & logFC < -0.07, "DE down p<0.05", "n.s.")))
limma_results$DE <- 
  factor(limma_results$DE, 
         ordered = TRUE, 
         levels = c("DE up p<0.05","DE down p<0.05","n.s."))

limma_results <- limma_results[order(limma_results$pvalue.limma,decreasing = F),]

customPlot2 <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","gray80")), 
  scale_colour_manual(values=c("#E43F3F","#268989","gray80"))
)

p3 <- ggplot(data = limma_results, aes(logFC, -log10(pvalue.limma), colour = DE)) +
  geom_vline(aes(xintercept = 0)) + xlim(-1,1)+ ylim(0,6.7)+
  geom_point(size=0.2) + geom_hline(aes(yintercept = -log10(0.05)),size=0.2) + geom_hline(aes(yintercept = -log10(1.17096e-05)),size=0.2,linetype="dashed")+
  geom_text_repel(aes(label = c(limma_results$Gene[limma_results$logFC>0][c(1:5)],limma_results$Gene[limma_results$logFC<0][c(1:5)]),segment.size=0.2), 
                  data = rbind(limma_results[limma_results$logFC>0,][c(1:5),],limma_results[limma_results$logFC<0,][c(1:5),]), 
                  vjust = 0, nudge_y = 0.1, size = 3)  +
  xlab("log2FoldChange") +ylab("")+
  customPlot+theme(legend.position = "none")+ggtitle("protein")+theme(plot.title = element_text(hjust = 0.5))

library(ggpubr)
p4<-ggarrange(p1,p2,p3,nrow=1,ncol=3)

ggsave("/path/to/Figures/1C_E.pdf",p4,width=12,height=4)

#### 1F+G #### - GO and KEGG enrichment analyses
library(fgsea)
library(msigdbr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)

# With DE genes p<0.05
prot <- limma_results$Gene[limma_results$pvalue.limma<0.05 & abs(limma_results$logFC)>0.07]
EIDs <- bitr(prot,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_KEGG <- EIDs[!duplicated(EIDs[c("ENTREZID")]),"ENTREZID"]

RNA <- DE_res$Gene[DE_res$pvalue<0.05 & abs(DE_res$log2FoldChange)>0.07]
EIDs <- bitr(RNA,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
RNA_KEGG <- EIDs[!duplicated(EIDs[c("ENTREZID")]),"ENTREZID"]

miR_targets_significant <- read.table("/path/to/miRNA/DE_downstream/miR_targets_significant.txt", quote="\"", comment.char="")
EIDs <- bitr(miR_targets_significant$V1,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
miR_target_KEGG <- EIDs[!duplicated(EIDs[c("ENTREZID")]),"ENTREZID"]

cc_GO <- compareCluster(geneClusters =list(miRNA=miR_targets_significant$V1,RNA=RNA,Protein=prot),fun = "enrichGO",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=0.05,ont="ALL")
cc_GO <- pairwise_termsim(cc_GO)

p4 <- emapplot(cc_GO,color="p.adjust",legend_n=2,cex_line=0.1,cex_label_category=0.8,layout="nicely",cex_category=0.7,showCategory=14)+scale_fill_manual(values=c("#bbbbbc","#0c70c8","#47a72f"))
ggsave("/path/to/Figures/1F.pdf",p4,height=5,width=6)

cc_KEGG <- compareCluster(geneClusters =list(miRNA=miR_target_KEGG,RNA=RNA_KEGG,Protein=prot_KEGG),fun = "enrichKEGG",keyType = "kegg",pvalueCutoff=0.05)
cc_KEGG <- pairwise_termsim(cc_KEGG)

p5 <- emapplot(cc_KEGG,color="p.adjust",legend_n=2,cex_line=0.1,cex_label_category=0.8,layout="nicely",cex_category=0.7,showCategory=18)+scale_fill_manual(values=c("#bbbbbc","#0c70c8","#47a72f"))
ggsave("/path/to/Figures/1G.pdf",p5,height=5,width=6)

##########################
#### Figure 2 - WGCNA #### 
##########################

#### 2A+B #### generated during network construction from net.rds object (see script 6b_WGCNA_automatic_network_construction.R)
# Example code for generating the dendrogram plot:
mergedColors = labels2colors(net$colors)
png("VS_dendogram.png", width = 6, height = 4,units="in",res=600)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Correlation heatmap created using moduleTraitCor, datTraits, and textMatrix (see script 6c_WGCNA_relate_to_traits.R)
# RNA CUD-associated modules
pdf("/path/to/Figures/2A_1.pdf", width = 9, height = 4)
par(mar = c(6, 8.5, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor[rownames(moduleTraitCor) %in% c("MEgrey60","MEskyblue3","MEmidnightblue","MElightcyan1"),c(1,6,2:5,7:12)],
               xLabels = names(datTraits)[c(1,6,2:5,7:12)],
               yLabels = c("MEgrey60","MEskyblue3","MEmidnightblue","MElightcyan1"),
               ySymbols = c("MEgrey60","MEskyblue3","MEmidnightblue","MElightcyan1"),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[c(12,19,21,40),c(1,6,2:5,7:12)],
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Protein CUD-associated modules 
pdf("/path/to/Figures/2A_2.pdf", width = 6.25, height = 3.5)
par(mar = c(6, 8.5, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor[rownames(moduleTraitCor) %in% c("MEyellow","MEbrown","MEtan"),c(1,6,2:5)],
               xLabels = names(datTraits)[c(1,6,2:5)],
               yLabels = c("MEyellow","MEbrown","MEtan"),
               ySymbols = c("MEyellow","MEbrown","MEtan"),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[c(6,7,16),c(1,6,2:5)],
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#### 2C ####
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# GO terms
# RNA module genes 
geneInfo0_Expr <- read.csv("/path/to/RNA/DE_downstream/WGCNA/geneInfo0_Expr.txt", sep=";")
# Protein module genes
geneInfo0 <- read.csv("/path/to/protein/data_analysis_results_V1/WGCNA/geneInfo_Prot0.txt", sep=";")

GO_gene_list2 <- list(yellow=geneInfo0$geneSymbol[geneInfo0$moduleColor == "yellow"],brown=geneInfo0$geneSymbol[geneInfo0$moduleColor == "brown"],tan=geneInfo0$geneSymbol[geneInfo0$moduleColor == "tan"],grey60=geneInfo0_Expr$geneSymbol[geneInfo0_Expr$moduleColor == "grey60"],skyblue3=geneInfo0_Expr$geneSymbol[geneInfo0_Expr$moduleColor == "skyblue3"],midnightblue=geneInfo0_Expr$geneSymbol[geneInfo0_Expr$moduleColor == "midnightblue"],lightcyan1=geneInfo0_Expr$geneSymbol[geneInfo0_Expr$moduleColor == "lightcyan1"])

module_GO2 <- compareCluster(geneClusters = GO_gene_list2,fun = "enrichGO",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=0.05,ont="ALL")
module_GO2 <-pairwise_termsim(module_GO2)

png("/path/to/Figures/2C.png",units="in",res=900,width=16,height=10,family="Arial")
treeplot(module_GO2,nCluster=9,geneClusterPanel="dotplot",offset=rel(2.5),nwords=10,  offset.params = list(bar_tree = rel(1), tiplab = rel(1.1), extend = 0, hexpand = 0.1),group_color=magma(15)[c(2:10)],hclust_method="ward.D",showCategory=6)
dev.off()

#### 2D+E ####
# GeneOverlap for WGCNA module genes, DE genes and marker genes, also miRNA target genes 
# Import module genes 
geneInfo0_Expr <- read.csv("/path/to/RNA/DE_downstream/WGCNA/geneInfo0_Expr.txt", sep=";")
gy_RNA <- unique(geneInfo0_Expr$genx[geneInfo0_Expr$moduleColor=="grey60"])
sb_RNA <- unique(geneInfo0_Expr$genx[geneInfo0_Expr$moduleColor=="skyblue3"])
mnb_RNA <- unique(geneInfo0_Expr$genx[geneInfo0_Expr$moduleColor=="midnightblue"])
cy_RNA <- unique(geneInfo0_Expr$genx[geneInfo0_Expr$moduleColor=="lightcyan1"])

geneInfo_Prot0 <- read.csv("/path/to/protein/data_analysis_results_V1/WGCNA/geneInfo_Prot0.txt", sep=";")
yl_Prot <- unique(geneInfo_Prot0$genx[geneInfo_Prot0$moduleColor=="yellow"])
br_Prot <- unique(geneInfo_Prot0$genx[geneInfo_Prot0$moduleColor=="brown"])
tan_Prot <- unique(geneInfo_Prot0$genx[geneInfo_Prot0$moduleColor=="tan"])

# Import DE genes and proteins
DE_RNA <- read.csv("/path/to/RNA/DE_results/DE_results.txt")
RNA_up <- unique(DE_RNA$Gene[DE_RNA$pvalue<0.05 & DE_RNA$log2FoldChange>0.07])
RNA_down <- unique(DE_RNA$Gene[DE_RNA$pvalue<0.05 & DE_RNA$log2FoldChange<0.07])

DE_Prot <- read.csv("/path/to/protein/data_analysis_results_V1/Limma_results_V1.csv")
Prot_up <- unique(DE_Prot$Gene[DE_Prot$pvalue.limma<0.05 & DE_Prot$logFC>0.07])
Prot_down <- unique(DE_Prot$Gene[DE_Prot$pvalue.limma<0.05 & DE_Prot$logFC<0.07])

# Create GeneOverlap matrix
query <- list(grey60_RNA = gy_RNA,skyblue3_RNA=sb_RNA,midnightblue_RNA = mnb_RNA,lightcyan1_RNA=cy_RNA,yellow_protein=yl_Prot,brown_protein=br_Prot,tan_protein=tan_Prot)
ref <- list(DE_up_RNA = RNA_up,DE_down_RNA=RNA_down,DE_up_protein=as.character(Prot_up),DE_down_protein=as.character(Prot_down))

library(GeneOverlap)  
gom.obj <- newGOM(ref,query,genome.size=19659)

# Plot GeneOverlap Figure
pdf("/path/to/Figures/2D_E.pdf",height=8,width=11)
drawHeatmap(gom.obj,what= "Jaccard", grid.col="Reds",note.col="black",adj.p = T,cutoff = 0.05)
dev.off()

#########################
#### Figure 3 - MOFA #### 
#########################

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

#### 3A ####
p1 <- plot_variance_explained(VS)
ggsave("/path/to/Figures/3A.pdf",p1, width = 4, height = 5)

#### 3B ####
p2 <- correlate_factors_with_covariates(VS, 
                                        covariates = colnames(pheno_cc[c("RIN","CUD","pH","Age","PMI","batch","Cocaine at death","OPC","Astrocyte","Oligodendrocyte","Inh_MSN" ,"Microglia","Inh_GABA")]), 
                                        plot="log_pval"
)
ggsave("/path/to/Figures/3B.pdf",p2, width = 3.5, height = 3.8)

#### 3C ####
p3 <- plot_factor(VS, 
                  factors = 10, 
                  color_by = "CUD",
                  add_violin = F,
                  dodge = TRUE
)+scale_fill_manual(values=c("#268989","#E43F3F"))+ylim(-2,3)
ggsave("/path/to/Figures/3C.pdf",p3, width = 4, height = 3)

#### 3D ####
library(cowplot)
p4 <- plot_top_weights(VS,
                       view = "miRNA",
                       factor = 10,
                       nfeatures = 10,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
)+panel_border(remove=T)

ggsave("/path/to/Figures/3D.pdf",p4, width = 4.2, height =3.5)

#### 3E ####
p5 <- plot_top_weights(VS,
                       view = "mRNA",
                       factor = 10,
                       nfeatures = 10,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
)+panel_border(remove=T)

ggsave("/path/to/Figures/3E.pdf",p5, width = 4.2, height = 3.5)

#### 3F ####
p6 <- plot_top_weights(VS,
                       view = "protein",
                       factor = 10,
                       nfeatures = 10,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
)+panel_border(remove=T)
ggsave("/path/to/Figures/3F.pdf",p6, width = 4.2, height = 3.5)

#### 3G ####
f10_mRNA <- get_weights(VS,views="mRNA")
f10_mRNA<-f10_mRNA$mRNA
rownames(f10_mRNA) <- gsub("_mRNA","",rownames(f10_mRNA))
f10_mRNA <- f10_mRNA[,10]
f10_mRNA <- sort(f10_mRNA,decreasing = TRUE)

f10_protein <- get_weights(VS,views="protein")
f10_protein<-f10_protein$protein
rownames(f10_protein) <- gsub("_protein","",rownames(f10_protein))
f10_protein <- f10_protein[,10]
f10_protein <- sort(f10_protein,decreasing = TRUE)

gene_list <- list(Factor10_mRNA=f10_mRNA ,Factor10_protein=f10_protein)

# Plot and write output

GO <- compareCluster(geneClusters = gene_list,fun = "gseGO",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=0.05,ont="ALL",exponent = 1,
                     minGSSize = 10,maxGSSize=500,eps=0,by="fgsea")

GO2 <- GO
GO2@compareClusterResult$Cluster <- as.character(GO2@compareClusterResult$Cluster)
GO2@compareClusterResult$Cluster[GO2@compareClusterResult$NES > 0] <- paste0(GO2@compareClusterResult$Cluster[GO2@compareClusterResult$NES > 0],"_positive")
GO2@compareClusterResult$Cluster[GO2@compareClusterResult$NES < 0] <-  paste0(GO2@compareClusterResult$Cluster[GO2@compareClusterResult$NES < 0],"_negative")
GO2@compareClusterResult$Cluster <- factor(GO2@compareClusterResult$Cluster,levels=c("Factor10_mRNA_positive","Factor10_protein_positive","Factor10_mRNA_negative","Factor10_protein_negative"))
GO2<-pairwise_termsim(GO2)
p1 <- emapplot(GO2,legend_n=3,nCluster=4,cex_line=0.1,cex_label_category=1.2,layout='nicely',cex_category=0.65,showCategory=10,pie="equal")+scale_fill_manual(values=c("#0c70c8","#8FC7F8","#B0E5A3","#47a72f"))
ggsave(paste0("/path/to/Figures/3G.pdf"),p1,height=9,width=10)

################################################
#### Figure 4 - snRNAseq DE and GO analyses ####
################################################

#### 4A ####  
seurat <- readRDS("/path/to/snRNAseq/2_data_processed/2_harmony_object.rds")
Idents(seurat)<-"celltypes"
p1 <- DimPlot(seurat, label=F, label.size=3,reduction = "umap.harmony",pt.size=0.000005,order=c("Lymphocyte","Endothelial","Ependymal","GABAergic-3","GABAergic-2","GABAergic-1","D2-MSN","D1-MSN","OPC","Oligodendrocyte","Microglia", "Astrocyte"),cols= c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#85CEE4","#006960","#003E65",  "#BFE1D7","#D1BCDC","#FCD8C1","#969997"))+NoAxes()
ggsave("/path/to/Figures/4A.pdf",plot = p1,width = 10,height=8)

#### 4B ####
genes <- c("GAD1","GAD2","DRD1","TAC1","DRD2","ADORA2A","HTR7","SST","PTHLH","PVALB","CCK","CALB2","AQP4","SLC1A2","VCAN","PDGFRA","MOBP","PLP1","CSF1R","P2RY12","CD96","CD3D","FLT1","CLDN5","CFAP299")
seurat@active.ident <- factor(seurat@active.ident,
                              levels= rev(c("D1-MSN",
                                            "D2-MSN",
                                            "GABAergic-1",
                                            "GABAergic-2",
                                            "GABAergic-3",
                                            "Astrocyte",
                                            "OPC",
                                            "Oligodendrocyte",
                                            "Microglia",
                                            "Lymphocyte",
                                            "Endothelial",
                                            "Ependymal")))

dp <- DotPlot(object = seurat, features = genes,assay ="RNA", dot.scale = 4.5) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=10)) + 
  theme(text=element_text(size = 8))+
  theme(aspect.ratio =0.8) + xlab("")+ylab("")
ggsave("/path/to/Figures/4B.pdf",plot = dp,width = 7,height=5)

#### 4C ####
# Import DE results - log2FC 0.5 cutoff and padj 0.05, only celltypes with more than 100 cells from each condition
OPC <- read.csv("/path/to/snRNAseq/3_results/DE/clusterOPC.csv")
OPC_sig <- OPC[OPC$p_val_adj<0.05 & abs(OPC$avg_log2FC)>0.25 ,]
OPC_sig <- OPC_sig[OPC_sig$pct.1 > 0.25 | OPC_sig$pct.2 > 0.25, ]
OPC_top10 <- OPC_sig$X[order(abs(OPC_sig$avg_log2FC),decreasing = T)][c(1:10)]

Oligodendrocyte <- read.csv("/path/to/snRNAseq/3_results/DE/clusterOligodendrocyte.csv")
Oligodendrocyte_sig <- Oligodendrocyte[Oligodendrocyte$p_val_adj<0.05 & abs(Oligodendrocyte$avg_log2FC)>0.25,]
Oligodendrocyte_sig <- Oligodendrocyte_sig[Oligodendrocyte_sig$pct.1 > 0.25 | Oligodendrocyte_sig$pct.2 > 0.25, ]
Oligodendrocyte_top10 <- Oligodendrocyte_sig$X[order(abs(Oligodendrocyte_sig$avg_log2FC),decreasing = T)][c(1:10)]

Astrocyte <- read.csv("/path/to/snRNAseq/3_results/DE/clusterAstrocyte.csv")
Astrocyte_sig <- Astrocyte[Astrocyte$p_val_adj<0.05 & abs(Astrocyte$avg_log2FC)>0.25,]
Astrocyte_sig <- Astrocyte_sig[Astrocyte_sig$pct.1 > 0.25 | Astrocyte_sig$pct.2 > 0.25, ]
Astrocyte_top10 <- Astrocyte_sig$X[order(abs(Astrocyte_sig$avg_log2FC),decreasing = T)][c(1:10)]

Microglia <- read.csv("/path/to/snRNAseq/3_results/DE/clusterMicroglia.csv")
Microglia_sig <- Microglia[Microglia$p_val_adj<0.05 & abs(Microglia$avg_log2FC)>0.25,]
Microglia_sig <- Microglia_sig[Microglia_sig$pct.1 > 0.25 | Microglia_sig$pct.2 > 0.25, ]
Microglia_top10 <- Microglia_sig$X[order(abs(Microglia_sig$avg_log2FC),decreasing = T)][c(1:10)]

GABAergic_1 <- read.csv("/path/to/snRNAseq/3_results/DE/clusterGABAergic-1.csv")
GABAergic_1_sig <- GABAergic_1[GABAergic_1$p_val_adj<0.05 & abs(GABAergic_1$avg_log2FC)>0.25,]
GABAergic_1_sig <- GABAergic_1_sig[GABAergic_1_sig$pct.1 > 0.25 | GABAergic_1_sig$pct.2 > 0.25, ]
GABAergic_1_top10 <- GABAergic_1_sig$X[order(abs(GABAergic_1_sig$avg_log2FC),decreasing = T)][c(1:10)]

D2MSN <- read.csv("/path/to/snRNAseq/3_results/DE/clusterD2-MSN.csv")
D2MSN_sig <- D2MSN[D2MSN$p_val_adj<0.05 & abs(D2MSN$avg_log2FC)>0.25,]
D2MSN_sig <- D2MSN_sig[D2MSN_sig$pct.1 > 0.25 | D2MSN_sig$pct.2 > 0.25, ]
D2MSN_top10 <- D2MSN_sig$X[order(abs(D2MSN_sig$avg_log2FC),decreasing = T)][c(1:10)]

D1MSN <- read.csv("/path/to/snRNAseq/3_results/DE/clusterD1-MSN.csv")
D1MSN_sig <- D1MSN[D1MSN$p_val_adj<0.05 & abs(D1MSN$avg_log2FC)>0.25,]
D1MSN_sig <- D1MSN_sig[D1MSN_sig$pct.1 > 0.25 | D1MSN_sig$pct.2 > 0.25, ]
D1MSN_top10 <- D1MSN_sig$X[order(abs(D1MSN_sig$avg_log2FC),decreasing = T)][c(1:10)]

# Get list of DE genes that need to be included in the heatmap
DE_genes <- c(OPC_sig$X,Oligodendrocyte_sig$X,D1MSN_sig$X,D2MSN_sig$X,GABAergic_1_sig$X,Microglia_sig$X,Astrocyte_sig$X)
DE_genes <- DE_genes[!duplicated(DE_genes)]

# Obtain log2FCs for each cluster for the selection of DE genes
OPC <- OPC[OPC$X %in% DE_genes,c("X","avg_log2FC")]
colnames(OPC)[2] <- "OPC"

Oligodendrocyte <- Oligodendrocyte[Oligodendrocyte$X %in% DE_genes,c("X","avg_log2FC")]
colnames(Oligodendrocyte)[2] <- "Oligodendrocyte"

Astrocyte <- Astrocyte[Astrocyte$X %in% DE_genes,c("X","avg_log2FC")]
colnames(Astrocyte)[2] <- "Astrocyte"

Microglia <- Microglia[Microglia$X %in% DE_genes,c("X","avg_log2FC")]
colnames(Microglia)[2] <- "Microglia"

D1MSN <- D1MSN[D1MSN$X %in% DE_genes,c("X","avg_log2FC")]
colnames(D1MSN)[2] <- "D1MSN"

D2MSN <- D2MSN[D2MSN$X %in% DE_genes,c("X","avg_log2FC")]
colnames(D2MSN)[2] <- "D2MSN"

GABAergic_1 <- GABAergic_1[GABAergic_1$X %in% DE_genes,c("X","avg_log2FC")]
colnames(GABAergic_1)[2] <- "GABAergic_1"

# merge dfs
library(tidyverse)
dfs<-list(OPC,Oligodendrocyte,Astrocyte,Microglia,D1MSN,D2MSN,GABAergic_1) %>% purrr::reduce(inner_join,by="X")
rownames(dfs)<-dfs$X

# Create log2FC heatmap
hm_DEG <- as.matrix(dfs[,-1])

library(ComplexHeatmap)
#reorder the rows to match with annotation labels
rows_DEG <- row_order(Heatmap(hm_DEG,name="log2FC",column_names_gp  = gpar(fontsize=8),column_names_rot = 45))
hm_DEG <- hm_DEG[rows_DEG,]

# subset genes for labeling in the heatmap
lab <- rownames(hm_DEG)
lab_genes <- c(OPC_top10,Oligodendrocyte_top10,Astrocyte_top10,Microglia_top10,D1MSN_top10,D2MSN_top10,GABAergic_1_top10)

#remove mitochondrial genes, pseudogenes and XIST
lab_genes <- lab_genes[-grep("RPS|RPL|XIST|MT-|AC|AL|MTRNR",lab_genes)]
lab[!(lab %in% lab_genes)] <- ""
lab_genes <- unique(lab_genes[!is.na(lab_genes)])

# Plot heatmap
library(ComplexHeatmap)
pdf("/path/to/Figures/4C.pdf",width=4,height=8)
Heatmap(hm_DEG,name="log2FC",cluster_columns = T,column_names_gp  = gpar(fontsize=8),column_names_rot = 45) + rowAnnotation(link = anno_mark(at = which(rownames(hm_DEG)==lab),                                                                                                                                     labels = lab[lab !=""], padding = unit(1, "mm"),link_gp = gpar(lwd = 0.1) ,labels_gp = gpar(fontsize = 6)))
dev.off()

#### 4D ####
# use the compareCluster function of clusterProfiler
# create gene list with DE genes for each cluster for up- and downregulated genes together
GO_both <- list(Astrocyte=unique(Astrocyte_sig$X),Microglia=unique(Microglia_sig$X),Oligodendrocyte=unique(Oligodendrocyte_sig$X),OPC=unique(OPC_sig$X),D1MSN=unique(D1MSN_sig$X),D2MSN=unique(D2MSN_sig$X),GABAergic_1=unique(GABAergic_1_sig$X))
both_GO <- compareCluster(geneClusters = GO_both,fun = "enrichGO",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=0.05)
both_GO <-pairwise_termsim(both_GO)
p <- emapplot(both_GO,legend_n=3,cex_line=0.1,cex_label_category=1.2,layout="nicely",cex_category=0.8,showCategory=10,pie="equal",repel=T)+scale_fill_manual(values=c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#85CEE4","#006960"))
ggsave("/path/to/Figures/4D.pdf",p,height=9,width=9)

#### 4E ####
OPC <- read.csv("/path/to/snRNAseq/3_results/DE/clusterOPC.csv")
Oligodendrocyte <- read.csv("/path/to/VS_multiome/snRNAseq/3_results/DE/clusterOligodendrocyte.csv")
Astrocyte <- read.csv("/path/to/snRNAseq/3_results/DE/clusterAstrocyte.csv")
Microglia <- read.csv("/path/to/snRNAseq/3_results/DE/clusterMicroglia.csv")
GABAergic_1 <- read.csv("/path/to/snRNAseq/3_results/DE/clusterGABAergic-1.csv")
D2MSN <- read.csv("/path/to/snRNAseq/3_results/DE/clusterD2-MSN.csv")
D1MSN <- read.csv("/path/to/snRNAseq/3_results/DE/clusterD1-MSN.csv")

# Obtain log2FCs and padj for each cluster 
OPC <- OPC[,c("X","avg_log2FC","p_val_adj")]
colnames(OPC)[c(2,3)] <- paste0(colnames(OPC)[c(2,3)],"_OPC")

Oligodendrocyte <- Oligodendrocyte[,c("X","avg_log2FC","p_val_adj")]
colnames(Oligodendrocyte)[c(2,3)] <- paste0(colnames(Oligodendrocyte)[c(2,3)],"_Oligodendrocyte")

Astrocyte <- Astrocyte[,c("X","avg_log2FC","p_val_adj")]
colnames(Astrocyte)[c(2,3)] <- paste0(colnames(Astrocyte)[c(2,3)],"_Astrocyte")

Microglia <- Microglia[,c("X","avg_log2FC","p_val_adj")]
colnames(Microglia)[c(2,3)] <- paste0(colnames(Microglia)[c(2,3)],"_Microglia")

D1MSN <- D1MSN[,c("X","avg_log2FC","p_val_adj")]
colnames(D1MSN)[c(2,3)] <- paste0(colnames(D1MSN)[c(2,3)],"_D1MSN")

D2MSN <- D2MSN[,c("X","avg_log2FC","p_val_adj")]
colnames(D2MSN)[c(2,3)] <- paste0(colnames(D2MSN)[c(2,3)],"_D2MSN")

GABAergic_1 <- GABAergic_1[,c("X","avg_log2FC","p_val_adj")]
colnames(GABAergic_1)[c(2,3)] <- paste0(colnames(GABAergic_1)[c(2,3)],"_GABAergic_1")

# merge dfs
library(tidyverse)
dfs<-list(Microglia,OPC,Oligodendrocyte,Astrocyte,D1MSN,D2MSN,GABAergic_1) %>% purrr::reduce(inner_join,by="X")
dfs_l2FC <- dfs[,c(1,2,4,6,8,10,12,14)]
dfs_p <- dfs[,c(1,3,5,7,9,11,13,15)]
colnames(dfs_l2FC) <- gsub("avg_log2FC_","",colnames(dfs_l2FC))
rownames(dfs_l2FC) <- dfs_l2FC$X
colnames(dfs_p) <- gsub("p_val_adj_","",colnames(dfs_p))
rownames(dfs_p) <- dfs_p$X
dfs_l2FC <- dfs_l2FC[,-1]
dfs_p <- dfs_p[,-1]

# structural component of ribosome genes
ribo <- gsub("/",";",both_GO@compareClusterResult[1,"geneID"])
ribo <- unlist(strsplit(ribo, ";"))

# electron transfer activity genes
nadh <- gsub("/",";",both_GO@compareClusterResult[3,"geneID"])
nadh <- unlist(strsplit(nadh, ";"))

# Create log2FC hms
hm_ribo_genes <- dfs_l2FC[rownames(dfs_l2FC) %in% ribo,]
hm_ribo_genes <- as.matrix(hm_ribo_genes[match(ribo,rownames(hm_ribo_genes)),])

hm_nadh_genes <- dfs_l2FC[rownames(dfs_l2FC) %in% nadh,]
hm_nadh_genes <- as.matrix(hm_nadh_genes[match(nadh,rownames(hm_nadh_genes)),])


## Create pval hms
makeStars <- function(x){
  stars <- c( "***", "**", "*","")
  vec <- c(0, 0.001, 0.01, 0.05,1.01)
  i <- findInterval(x, vec)
  stars[i]
}

hm_ribo_genes_p <- dfs_p[rownames(dfs_p) %in% ribo,]
hm_ribo_genes_p <- as.matrix(hm_ribo_genes_p[match(ribo,rownames(hm_ribo_genes_p)),])
hm_ribo_genes_star <-  t(apply(hm_ribo_genes_p,1,makeStars))
colnames(hm_ribo_genes_star) <-  colnames(hm_ribo_genes_p)

hm_nadh_genes_p <- dfs_p[rownames(dfs_p) %in% nadh,]
hm_nadh_genes_p <- as.matrix(hm_nadh_genes_p[match(nadh,rownames(hm_nadh_genes_p)),])
hm_nadh_genes_star <-  t(apply(hm_nadh_genes_p,1,makeStars))
colnames(hm_nadh_genes_star) <-  colnames(hm_nadh_genes_p)

library(RColorBrewer)
library(pheatmap)

pdf("/path/to/Figures/4E_1.pdf",width=6,height=12)
pheatmap(hm_ribo_genes, breaks=seq(-1.5, 1.5, length.out=10),color = colorRampPalette(c("blue","white","red"))(10),cluster_rows = F,cluster_cols = F,fontsize = 12,display_numbers = hm_ribo_genes_star,number_color="black",angle_col = 45,cellwidth = 20,cellheight = 20,main = "ribosomal genes")
dev.off()

# ETC includes genes from NADH dehydrogenase and cytochrome c oxidase 
pdf("/path/to/Figures/4E_2.pdf",width=6,height=12)
pheatmap(hm_nadh_genes, breaks=seq(-1.5, 1.5, length.out=10),color = colorRampPalette(c("blue","white","red"))(10),cluster_rows = F,cluster_cols = F,fontsize = 12,display_numbers = hm_nadh_genes_star,number_color="black",angle_col = 45,cellwidth = 20,cellheight = 20,main = "ETC genes")
dev.off()

######################################
#### Figure 5 - consensus hdWGCNA ####
######################################

#### 5A ####
seurat_merged <- readRDS("/path/to/Phillips_VS_data/2_data_processed/4_consensus_human_rat.rds")
# Plot UMAP of datasets
Idents(seurat_merged) <- "Species"
seurat_merged_rat <- subset(seurat_merged,idents=c("rat"))
seurat_merged_human <- subset(seurat_merged,idents=c("human"))

p0 <- DimPlot(seurat_merged, group.by='celltypes', raster=FALSE, pt.size=0.000005,order=c("Inh_GABA","Inh_MSN","OPC","Oligodendrocyte","Microglia", "Astrocyte"),cols= c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#006960"),reduction = "umap.harmony") + umap_theme()+ggtitle("Multi-species integration")
p1 <- DimPlot(seurat_merged_rat, group.by='celltypes', raster=FALSE, pt.size=0.000005,order=c("Inh_GABA","Inh_MSN","OPC","Oligodendrocyte","Microglia", "Astrocyte"),cols= c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#006960"),reduction = "umap.harmony") + umap_theme()+ggtitle("Rat - Phillips et al.")
p2 <- DimPlot(seurat_merged_human, group.by='celltypes', raster=FALSE, pt.size=0.000005,order=c("Inh_GABA","Inh_MSN","OPC","Oligodendrocyte","Microglia", "Astrocyte"),cols= c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#006960"),reduction = "umap.harmony") + umap_theme()+ggtitle("Human - DBCBB")
p3 <- p0+p2+p1

pdf("/path/to/Figures/5A.pdf",height=5, width=17)
p3
dev.off()

#### 5B ####
p_Inh_MSN <- ModuleFeaturePlot(seurat_merged, order='shuffle', raster=TRUE, raster_dpi=100, alpha=1,raster_scale=0.25, wgcna_name='Inh_MSN_consensus', restrict_range=FALSE,reduction = "umap.harmony")

pdf("/path/to/Figures/5B_1.pdf",height=20, width=20)
wrap_plots(c(list(p1),list(p2),list(p3),p_Inh_MSN), ncol=5)
dev.off()

p_Astrocyte <- ModuleFeaturePlot(seurat_merged, order='shuffle', raster=TRUE, raster_dpi=100, alpha=1, raster_scale=0.25, wgcna_name='Astrocyte_consensus', restrict_range=FALSE,reduction = "umap.harmony")

pdf("/path/to/Figures/5B_2.pdf",height=20, width=20)
wrap_plots(c(list(p1),list(p2),list(p3),p_Astrocyte), ncol=5)
dev.off()

#### 5C ####
seurat_human <- readRDS("/path/to/snRNAseq/2_data_processed/3_hdWGCNA_object.rds")
# Astrocyte
ast_merged <- GetModules(seurat_merged,wgcna_name="Astrocyte_consensus") 
ast_human <- GetModules(seurat_human,wgcna_name="Astrocyte") 

ast_list_merged <-list(Astrocyte_CM4=ast_merged$gene_name[ast_merged$module=="Astrocyte-CM4"],Astrocyte_CM6=ast_merged$gene_name[ast_merged$module=="Astrocyte-CM6"],Astrocyte_CM8=ast_merged$gene_name[ast_merged$module=="Astrocyte-CM8"])
ast_list_human <- list(Astrocyte_M12=ast_human$gene_name[ast_human$module=="Astrocyte-M12"],Astrocyte_M14=ast_human$gene_name[ast_human$module=="Astrocyte-M14"],Astrocyte_M18=ast_human$gene_name[ast_human$module=="Astrocyte-M18"])

library(GeneOverlap)
gom.obj <- newGOM(ast_list_human,ast_list_merged ,36588)

pdf("/path/to/Figures/5C.pdf",height=8, width=8)
drawHeatmap(gom.obj,grid.col="Reds", note.col="black",adj.p = T,what=c("Jaccard"))
dev.off()

#### 5D ####
# Inh_MSN
MSN_merged <- GetModules(seurat_merged,wgcna_name="Inh_MSN_consensus") 
MSN_human <- GetModules(seurat_human,wgcna_name="Inh_MSN") 

MSN_list_merged <-list(Inh_MSN_CM1=MSN_merged$gene_name[MSN_merged$module=="Inh_MSN-CM1"],Inh_MSN_CM4=MSN_merged$gene_name[MSN_merged$module=="Inh_MSN-CM4"],Inh_MSN_CM5=MSN_merged$gene_name[MSN_merged$module=="Inh_MSN-CM5"],Inh_MSN_CM7=MSN_merged$gene_name[MSN_merged$module=="Inh_MSN-CM7"],Inh_MSN_CM13=MSN_merged$gene_name[MSN_merged$module=="Inh_MSN-CM13"])
MSN_list_human <- list(Inh_MSN_M1=MSN_human$gene_name[MSN_human$module=="Inh_MSN-M1"],Inh_MSN_2=MSN_human$gene_name[MSN_human$module=="Inh_MSN-M2"],Inh_MSN_M7=MSN_human$gene_name[MSN_human$module=="Inh_MSN-M7"],Inh_MSN_M10=MSN_human$gene_name[MSN_human$module=="Inh_MSN-M10"])

gom.obj <- newGOM( MSN_list_human,MSN_list_merged ,36588)

pdf("/path/to/Figures/5D.pdf",height=8, width=9.5)
drawHeatmap(gom.obj,grid.col="Reds", note.col="black",adj.p = T,what=c("Jaccard"))
dev.off()

#### 5E ####
# Make GO/KEGG enrichment plots 
dbs  <- c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021","KEGG_2021_Human")
seurat_merged <- RunEnrichr(seurat_merged, dbs=dbs, wgcna_name="Astrocyte_consensus")
seurat_merged <- RunEnrichr(seurat_merged, dbs=dbs, wgcna_name="Inh_MSN_consensus")

# Astrocyte
terms_ast<-GetEnrichrTable(seurat_merged,wgcna_name = "Astrocyte_consensus")
terms_ast_GO<-terms_ast[terms_ast$db %in% c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021"),]
terms_ast_KEGG<-terms_ast[terms_ast$db %in% c("KEGG_2021_Human"),]

#  Astrocyte-CM4
m4a <- terms_ast_GO[terms_ast_GO$module=="Astrocyte-CM4",]
m4a_full <- m4a[order(m4a$Combined.Score,decreasing = T),][c(1:10),]
m4a <- m4a[order(m4a$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m4ak <- terms_ast_KEGG[terms_ast_KEGG$module=="Astrocyte-CM4",]
m4_ak_full <- m4ak[order(m4ak$Combined.Score,decreasing = T),][c(1:10),]
m4ak <- m4ak[order(m4ak$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m4a <- rbind(m4a,m4ak)
m4a$Term[c(1:10)] <- gsub(".{13}$","",m4a$Term[c(1:10)])
m4a$db[c(1:10)]<-"GO"
m4a$db[c(11:20)]<-"KEGG"

m4a$Combined.Score <- log(m4a$Combined.Score)
position1 <- rev(m4a$Term)

p1 <- ggplot(data = m4a[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="brown") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(3,8)
p2 <- ggplot(data = m4a[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="brown") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(3,7)

pdf("/path/to/Figures/5E.pdf",width=14,height=9)
plot_grid(p1,p2,ncol=1,align="v")
dev.off()

#### 5F ####
# Astrocyte-CM6
m6 <- terms_ast_GO[terms_ast_GO$module=="Astrocyte-CM6",]
m6_full <- m6[order(m6$Combined.Score,decreasing = T),][c(1:10),]
m6 <- m6[order(m6$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m6k <- terms_ast_KEGG[terms_ast_KEGG$module=="Astrocyte-CM6",]
m6_k_full <- m6k[order(m6k$Combined.Score,decreasing = T),][c(1:10),]
m6k <- m6k[order(m6k$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m6 <- rbind(m6,m6k)
m6$Term[c(1:10)] <- gsub(".{13}$","",m6$Term[c(1:10)])
m6$db[c(1:10)]<-"GO"
m6$db[c(11:20)]<-"KEGG"
m6$Combined.Score <- log(m6$Combined.Score)
position1 <- rev(m6$Term)

p1 <- ggplot(data = m6[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="red") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(3,8)
p2 <- ggplot(data = m6[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="red") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(3,7)

pdf("/path/to/Figures/5F.pdf",width=10,height=9)
plot_grid(p1,p2,ncol=1,align="v")
dev.off()

#### 5G ####
# Astrocyte-CM8
m8a <- terms_ast_GO[terms_ast_GO$module=="Astrocyte-CM8",]
m8a_full <- m8a[order(m8a$Combined.Score,decreasing = T),][c(1:10),]
m8a <- m8a[order(m8a$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m8ak <- terms_ast_KEGG[terms_ast_KEGG$module=="Astrocyte-CM8",]
m8_ak_full <- m8ak[order(m8ak$Combined.Score,decreasing = T),][c(1:10),]
m8ak <- m8ak[order(m8ak$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m8a <- rbind(m8a,m8ak)
m8a$Term[c(1:10)] <- gsub(".{13}$","",m8a$Term[c(1:10)])
m8a$db[c(1:10)]<-"GO"
m8a$db[c(11:20)]<-"KEGG"
m8a$Combined.Score <- log(m8a$Combined.Score)
position1 <- rev(m8a$Term)

p1 <- ggplot(data = m8a[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="black") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(3,8)
p2 <- ggplot(data = m8a[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="black") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(3,7)

pdf("/path/to/Figures/5G.pdf",width=10,height=9)
plot_grid(p1,p2,ncol=1,align="v")
dev.off()

#### 5H ####
# Inh_MSN-CM4
m4 <- terms_msn_GO[terms_msn_GO$module=="Inh_MSN-CM4",]
m4_full <- m4[order(m4$Combined.Score,decreasing = T),][c(1:10),]
m4 <- m4[order(m4$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m4k <- terms_msn_KEGG[terms_msn_KEGG$module=="Inh_MSN-CM4",]
m4_k_full <- m4k[order(m4k$Combined.Score,decreasing = T),][c(1:10),]
m4k <- m4k[order(m4k$Combined.Score,decreasing = T),][c(1:10),c(1,8,10)]
m4 <- rbind(m4,m4k)
m4$Term[c(1:10)] <- gsub(".{13}$","",m4$Term[c(1:10)])
m4$db[c(1:10)]<-"GO"
m4$db[c(11:20)]<-"KEGG"
m4$Combined.Score <- log(m4$Combined.Score)
position1 <- rev(m4$Term)

p1 <- ggplot(data = m4[c(1:10),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="blue") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(11:20)]) + theme(text=element_text(size=20))+xlab("")+xlim(3,8)
p2 <- ggplot(data = m4[c(2:20),], aes(x = Combined.Score, y = Term)) + geom_point(size=3,color="blue") +theme_bw() + ylab("") + xlab("Enrichment log(combined score)") + scale_y_discrete(limits = position1[c(1:10)]) + theme(text=element_text(size=20))+xlab("Enrichment log(combined score)")+xlim(3,7)

pdf("/path/to/Figures/5H.pdf",width=12,height=9)
plot_grid(p1,p2,ncol=1,align="v")
dev.off()

#############################
#### Figure 6 - CellChat ####
#############################

# For generation of CellChat objects see script 6_CellChat.Rmd
# Merge the two CellChat objects
cellchat <- readRDS("3_results/NeuronChat/CellChat_CUD.rda") # CUD
cellchat2 <-readRDS("3_results/NeuronChat/CellChat_Ctrl.rda") # Ctrl

object.list <- list(CUD = cellchat, Ctrl = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# DE analysis of ligand-receptor interactions

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "CUD"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in CUD
net.up <- subsetCommunication(cellchat, net = net, datasets = "CUD",ligand.logFC = 0.07, receptor.logFC = 0.00001)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in Ctrl, i.e.,downregulated in CUD
net.down <- subsetCommunication(cellchat, net = net, datasets = "Ctrl",ligand.logFC = -0.07, receptor.logFC = -0.00001)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Chord diagrams for differentially expressed ligand-receptor pairs

#### 6A ####
pdf("/path/to/Figures/6A_1.pdf",width=10,height=10)
netVisual_chord_gene(object.list[[1]], slot.name = 'net', net = net.up, lab.cex = 1, big.gap = 100,small.gap = 100,color.use=c("#F7AB64","#3A9BCC","#85CEE4"),reduce=0.0025)
dev.off()

#### 6B ####
pdf("/path/to/Figures/6B_1.pdf",width=12,height=12)
netVisual_chord_gene(object.list[[2]], slot.name = 'net', net = net.down, lab.cex = 1, big.gap = 100,small.gap = 100,color.use=c("#F7AB64","#3A9BCC","#85CEE4"),reduce=0.0025)
dev.off()

# Pie plot of annotation categories 
up <- data.frame(Annotation=names(table(net.up$annotation)),Count=as.vector(unname(table(net.up$annotation))))
down <- data.frame(Annotation=names(table(net.down$annotation)),Count=as.vector(unname(table(net.down$annotation))))

up$perc <- up$Count/sum(up$Count)*100
up$perc2 <- paste0(as.character(round(up$perc,digits=1)),"%")
down$perc <- down$Count/sum(down$Count)*100
down$perc2 <- paste0(as.character(round(down$perc,digits=1)),"%")

# Add label position
up <- up %>%
  arrange(desc(Annotation)) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)

down <- down %>%
  arrange(desc(Annotation)) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)

p3 <- ggplot(up, aes(x = 2, y = perc, fill = Annotation)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = perc2), color = "white",size=5)+
  theme_void()+scale_fill_manual(values = magma(8)[c(2,3,5)]) +
  xlim(0.5, 2.5)


p4 <- ggplot(down, aes(x = 2, y = perc, fill = Annotation)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = perc2), color = "white",size=5)+
  theme_void()+scale_fill_manual(values = magma(8)[c(2,3,6,5)]) +
  xlim(0.5, 2.5)

ggsave("/path/to/Figures/6A_2.pdf",p3,width=6,height=6)
ggsave("/path/to/Figures/6B_2.pdf",p4,width=6,height=6)

