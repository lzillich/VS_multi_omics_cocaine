#### Correlation analysis of transcriptome and proteome ####
# Last modification 2024-12-04 EZ
# Analysis following the protocol from Yang & Gorski 2022, STAR Protocols (PMID: 35634361)

# Create transcript to gene matching based on the ENCODE reference transcriptome 
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file="/path/to/reference/gencode.v43.primary_assembly.annotation.gtf",format="auto",organism = "Homo sapiens")
k <- keys(txdb,keytype="TXNAME")
tx2gene <- select(txdb,k,"GENEID","TXNAME")

# Get salmon quant.sf files
setwd("/path/to/RNA/2_salmon_quant/")
salmon_files <- list.files(path="/path/to/RNA/2_salmon_quant/",pattern = "\\.sf$",recursive = T)
names(salmon_files) <- gsub("/quant.sf","",salmon_files)

# Import TPM using tximport
library(tximport)
txi <- tximport(salmon_files,type="salmon",tx2gene=tx2gene)
tpm <- txi$abundance

# Gene Annotation
required.packages <- c("tidyverse", "reshape2","biomaRt", "org.Hs.eg.db", "edgeR", "limma", "NMF", "sva", "vsn", "ggpubr", "ggfortify", "GSEABase")
sapply(c(required.packages, "magrittr"), require, character.only = T)

mart <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")

RNA.anno <- getBM(attributes = c("ensembl_gene_id_version", "hgnc_symbol", "entrezgene_id", "gene_biotype"), filters = "ensembl_gene_id_version", values = rownames(tpm), mart = mart, useCache = F) %>%
set_colnames(c("EnsemblGeneID", "HGNCSymbol", "EntrezID", "Biotype")) %>%
  
  dplyr::mutate(EntrezID = as.character(EntrezID)) %>%
  
  group_by(EnsemblGeneID, Biotype) %>%
  
  dplyr::summarise(HGNCSymbol = ifelse(length(unique(HGNCSymbol)) > 1, paste(unique(HGNCSymbol), collapse = ", "), unique(HGNCSymbol)), EntrezID = ifelse(length(unique(EntrezID)) > 1, paste(unique(EntrezID), collapse = ", "), unique(EntrezID))) %>%
  
  ungroup %>%
  
  mutate_all(function(x) na_if(x, "NA")) %>%
  
  mutate_all(function(x) na_if(x, ""))

tpm <- as.data.frame(tpm)
tpm$EnsemblGeneID <- rownames(tpm)
tpm <- left_join(tpm, RNA.anno, by = "EnsemblGeneID")
tpm <- tpm[,c(41:44,1:40)]

# Create log2 transformed TPM values 
logTPM <- tpm %>% filter(., Biotype == "protein_coding") %>%
  column_to_rownames("EnsemblGeneID") %>%
  dplyr::select(matches("DE")) %>%
  as.matrix %>% {
    mat = .
    mat %>%
      is_greater_than(1) %>%
      rowSums %>%
      is_weakly_greater_than(ncol(mat)*0.5) %>%
      mat[.,] %>% return
  } %>%
  add(1) %>%
  log2

# Add gene symbols
logTPM <- as.data.frame(logTPM)
logTPM$EnsemblGeneID <- rownames(logTPM)
logTPM <- merge(RNA.anno,logTPM,by="EnsemblGeneID")
logTPM <- logTPM[,c(3,5:44)]

# Convert Seq IDs to Brain IDs
pheno <- read.csv("/path/to/RNA/DE_results/pheno_coc_at_death.txt", sep=";")
logTPM <- logTPM[!is.na(logTPM$HGNCSymbol),]
logTPM <- logTPM[!duplicated(logTPM$HGNCSymbol),]
rownames(logTPM) <- logTPM$HGNCSymbol
logTPM <- logTPM[,-1]

phe_ordered <- pheno[match(colnames(logTPM),pheno$rn),]
colnames(logTPM) <- phe_ordered$Brain_ID

# Get proteomics data (vsn,normalized,log2 transformed and batch corrected as used for WGCNA/MOFA)
datProt <- read.csv("/path/to/protein/data_analysis_results_V1/normalized_data_batch_removed.txt", row.names=1, sep=";")
rownames(datProt)<-gsub("CASE_rep_","",rownames(datProt))
rownames(datProt)<-gsub("CTRL_rep_","",rownames(datProt))
rownames(datProt)<-gsub("_0","",rownames(datProt))
rownames(datProt)<-gsub("_1","",rownames(datProt))
datProt <- t(datProt)

# Get intersection of genes 
shared <- intersect(rownames(datProt),rownames(logTPM))
samp <- intersect(colnames(datProt),colnames(logTPM))
logTPM_int <- logTPM[rownames(logTPM) %in% shared,colnames(logTPM) %in% samp]
datProt_int <- as.data.frame(datProt[rownames(datProt) %in% shared,colnames(datProt) %in% samp])

# order by rownames and colnames
logTPM_int <- logTPM_int[match(rownames(datProt_int),rownames(logTPM_int)),]
logTPM_int <- logTPM_int[,match(colnames(datProt_int),colnames(logTPM_int))]

# test for identical names
identical(rownames(logTPM_int),rownames(datProt_int))#TRUE
identical(colnames(logTPM_int),colnames(datProt_int))#TRUE
write.table(logTPM_int,"/path/to/RNA/mRNA_prot_corr/logTPM.txt",sep=" ",quote=F,row.names = T,col.names = T)
write.table(datProt_int,"/path/to/RNA/mRNA_prot_corr/logTMTint.txt",sep=" ",quote=F,row.names = T,col.names = T)

# Calculate the correlation between expr and prot

diag <- diag(cor(logTPM_int, datProt_int, method = "pearson")) 
 
  df_cor <- data.frame(Sample=names(diag),R=unname(diag))
  
  p1 <- ggdotchart(df_cor, x = "Sample", y = "R", ggtheme = theme_pubr(base_size = 10), add = "segment", ylim = c(0, 1), ylab = "mRNA-protein correlation (Pearson R)") +theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+xlab("Samples (N=39)")

pheno_with_CUD_Ctrl <- read.csv("/path/to/RNA/DE_results/pheno_with_CUD_Ctrl.txt", sep=";")
df_cor2 <-  df_cor
df_cor2$Sample <- gsub("X","",df_cor2$Sample)
colnames(df_cor2)[1]<-"Brain_ID"  
df_cor2 <- merge(df_cor2,pheno_with_CUD_Ctrl,by="Brain_ID")
write.table(df_cor2,"/path/to/RNA/mRNA_prot_corr/sample_corr_table.txt",sep=" ",quote=F,row.names = T,col.names = T)

# Test significance of correlation coefficients between CUD and Ctrl 
  
  df_cor$CUD <- phe_ordered$CUD[match(df_cor$Sample,phe_ordered$Brain_ID)] 
  wilcox.test(df_cor$R[df_cor$CUD ==1],df_cor$R[df_cor$CUD ==0])#p-value = 0.478
  
  # Comparison at the gene level 
  diag2 <- diag(cor(t(logTPM_int), t(datProt_int), method = "pearson"))
  
  diag2_df <- data.frame(Gene=names(diag2),R=unname(diag2))
  diag2_df <- diag2_df[order(diag2_df$R,decreasing = T),]
  write.table(diag2_df,"/path/to/RNA/mRNA_prot_corr/gene_corr_table.txt",sep=" ",quote=F,row.names = T,col.names = T)
  
  # Add median
  med <- median(diag2)
 
 # add quantiles 1% and 99% for gene annotation 
 q <- quantile(diag2, probs = c(0.01,0.5,0.99))

# Add labels of genes at both ends of distribution
 low <- sort(diag2[diag2<unname(q[1])])
 high <- sort(diag2[diag2>unname(q[3])])
 write.table(low,"/path/to/RNA/mRNA_prot_corr/corr_low_genes.txt",sep=" ",quote=F,row.names = T,col.names = F)
 write.table(high,"/path/to/RNA/mRNA_prot_corr/corr_high_genes.txt",sep=" ",quote=F,row.names = T,col.names = F)
 
 p2 <- gghistogram(diag2, bins = 50, theme = theme_pubr(base_size = 10), fill = "gray",
                   xlab = "mRNA-protein correlation (Pearson R)",ylab = "Number of genes") +xlim(-1,1)+geom_vline(xintercept=med,col="red")+geom_vline(xintercept=unname(q[1]),col="red",linetype="dashed")+geom_vline(xintercept=unname(q[3]),col="red",linetype="dashed")
 
 p3 <- ggarrange(p1,p2,ncol=1,nrow=2)
 
 ggsave("/path/to/RNA/mRNA_prot_corr/Corr_Figure.pdf",p3,width=6,height=8)
 
# Overall correlation plot 
 RNA_mean <- rowMeans(logTPM_int)
 prot_mean <- rowMeans(datProt_int)
 dat <- data.frame(RNA=unname(RNA_mean),Protein=unname(prot_mean),row.names = names(RNA_mean))
 
 p4 <- ggscatter(dat, x = "RNA", y = "Protein", 
                 size=0.1,color="gray10",
           cor.coef = TRUE, fullrange = T,cor.method = "pearson",
           xlab = "log2 TPM - RNA", ylab = "log2 TMT intensity - Protein",title="Correlation by gene")+theme(plot.title = element_text(face = "bold"))
 ggsave("/path/to/RNA/mRNA_prot_corr/Corr_mean_data.pdf",p4,width=6,height=6)
 
 
 # Select highly expressed genes 
 quantile(dat$RNA,probs=c(0,0.9,0.95,0.99,1))
 #0%        90%        95%        99%       100% 
 #0.9420165  5.4077049  6.1732164  7.3782927 10.0294791
 
 quantile(dat$Protein,probs=c(0,0.9,0.95,0.99,1))
 #      0%      90%      95%      99%     100% 
 #13.22876 21.87608 22.78064 24.52323 26.75720 
 
 # Top 40 are top 1% of expressed proteins - show the correlation for the highly expressed proteins only
 dat_high_expr <- dat[order(dat$Protein,decreasing = T),]
 dat_high_expr <- dat_high_expr[c(1:40),]
 p4_1 <- ggscatter(dat_high_expr, x = "RNA", y = "Protein", 
                 size=0.1,color="gray10",label = rownames(dat_high_expr),font.label = c(8, "plain"),repel = T,
                 cor.coef = TRUE, cor.method = "spearman",
                 xlab = "log2 TPM - RNA", ylab = "log2 TMT intensity - Protein",title="Highly expressed proteins")+theme(plot.title = element_text(face = "bold"))
 ggsave("/path/to/RNA/mRNA_prot_corr/Corr_mean_data_highly_expressed.pdf",p4_1,width=6,height=6)
 
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
 
# Create final plots for the CUD/Ctrl correlation coefficient comparison
library(ggpubr)

p7 <- ggarrange(p4,p4_1,nrow=2,ncol=1)
ggsave("/path/to/RNA/mRNA_prot_corr/Corr_plots_combined.pdf",p7,width=6,height=12)

dat_CUD_comb$status <- "CUD"
rownames(dat_CUD_comb) <- paste0("CUD_",rownames(dat_CUD_comb))
dat_Ctrl_comb$status <- "Ctrl"
rownames(dat_Ctrl_comb) <- paste0("Ctrl_",rownames(dat_Ctrl_comb))

dat_combined <- rbind(dat_CUD_comb,dat_Ctrl_comb)

p7_1 <- ggscatter(dat_combined , x = "RNA", y = "Protein", 
                size=0.1,color = "status",
                palette = c("gray10","red"),
                cor.coef = F, fullrange = T,cor.method = "pearson",
                xlab = "log2 TPM - RNA", ylab = "log2 TMT intensity - Protein",title="Correlation by CUD status")+theme(plot.title = element_text(face = "bold"))+stat_cor(aes(color = status))+theme(legend.position = "none")
ggsave("/path/to/RNA/mRNA_prot_corr/Corr_plots_CUD_Ctrl.pdf",p7_1,width=6,height=6)


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

# GSEA on the correlation coefficients to infer pathways

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
gse_neg2 <- gse_neg2[c(1:10),]
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

ggsave("/path/to/RNA/mRNA_prot_corr/GSEA_corr_waterfall.pdf",gsea_plot,width=12,height=8)
