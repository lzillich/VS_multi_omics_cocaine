## Data Preparation WGCNA 
## Code based on WGCNA tutorials by Langfelder & Horvath
# Adapted by Lea Zillich and Eric Zillich 
# last change: EZ 2024-12-04

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

#mean(pheno$pH,na.rm=T) 6.335
pheno$pH[is.na(pheno$pH)] <- 6.335
pheno_cc <- merge(pheno,pheno_c,by="Brain_ID")

# remove PCA outliers 
pheno_cc <- pheno_cc[!(pheno_cc$rn %in% c("DE80NGSUKBR136473","DE42NGSUKBR136478","DE86NGSUKBR136462")),]

  # Load the expression and trait data saved in the first part
  setwd("/path/to/DE_downstream/WGCNA/")
  lnames = load(file ="VS_Expr_dataInput_processed_filtered.RData") 
  #The variable lnames contains the names of loaded variables.
  lnames

  # Load network data saved in the second part.
  lnames = load(file = "/path/to/miRNA/VS_networkConstruction_Expr.RData");
  lnames
  
  # Define numbers of genes and samples
  setwd(dir0)
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p",method = "pearson")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  write.table(MEs, file = "/path/to/DE_downstream/WGCNA/MEs_Expr_filtered.txt", sep = ";", quote = F, row.names = T)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  # Display the correlation values within a heatmap plot
  names(datTraits)[6] <- "Cocaine at death"
  
  pdf("/path/to/miRNA/DE_downstream/WGCNA/Heatmap_modules_miR_Expr_filtered.pdf", width = 6.5, height= 6.5)
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
  
  
  # Define variable weight containing the weight column of datTrait
  CUD = as.data.frame(datTraits$CUD)
  names(CUD) = "CUD"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, CUD, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(CUD), sep="");
  names(GSPvalue) = paste("p.GS.", names(CUD), sep="")
  
  mod <- DF(moduleTraitCor)
  module <- rownames(mod)[mod$CUD == max(abs(mod$CUD))]
  module <- substr(module, 3, nchar(module))
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  
  png("/path/to/miRNA/DE_downstream/WGCNA/module_membership_gene_significance_Expr.png",height=12,width=12,res=300,units="in")
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for CUD status",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 3, cex.lab = 3, cex.axis = 3, col = module)
  dev.off()
  
  write.table(mod, "/path/to/miRNA/DE_downstream/WGCNA/mod_trait_cor_Expr.txt", sep = ";", quote = F)

#### RNA data ####

dir0 <- "/path/to/RNA"
  
  setwd(dir0)
  
  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE);
  
  # Load the expression and trait data saved in the first part
  setwd("DE_downstream/WGCNA/")
  lnames = load(file ="VS_Expr_dataInput_processed_filtered.RData") 
  #The variable lnames contains the names of loaded variables.
  lnames
  
  # Load network data saved in the second part.
  lnames = load(file = "/path/to/RNA/VS_networkConstruction_Expr.RData");
  lnames
  
  # Define numbers of genes and samples
  setwd(dir0)
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  write.table(MEs, file = "DE_downstream/WGCNA/MEs_Expr_filtered.txt", sep = ";", quote = F, row.names = T)
  
  # Regression analysis for the modules showing significant correlation with CUD status
  ME_reg_dat <- merge(MEs,datTraits, by=0)
  ME_reg_dat$CUD <- as.factor(ME_reg_dat$CUD)
  mods <- colnames(MEs)
  beta_mat <- matrix(0, nrow=5,ncol=length(mods))
  colnames(beta_mat)<-mods
  rownames(beta_mat) <- colnames(datTraits)[1:5]
  p_mat <- matrix(0, nrow=5,ncol=length(mods))
  colnames(p_mat)<-mods
  rownames(p_mat) <- colnames(datTraits)[1:5]
  
  for(i in mods){
    lm <- lm(get(i)~CUD+Age+PMI+pH+RIN, data = ME_reg_dat)
    beta_mat[,i] <- summary(lm)$coefficients[-1,1]
    p_mat[,i] <-summary(lm)$coefficients[-1,4]
  }
  
  p_mat[1,][p_mat[1,] <0.05 ] # modules significantly associated with CUD status in the regression analysis
  #MEivory        MEgrey60 MEpaleturquoise      MEskyblue3  MEmidnightblue         MEbrown    MElightcyan1        MEpurple 
  #0.011115386     0.002098758     0.027848830     0.004279865     0.015366851     0.024717044     0.036322610     0.020744202
  
  colnames(p_mat)<-paste0("p_",colnames(p_mat))
  reg_res <- as.data.frame(matrix(0,nrow=5,ncol=108))
  colnames(reg_res)<- c(rbind(colnames(beta_mat),colnames(p_mat)))
  rownames(reg_res)<-rownames(beta_mat)
  
  
  for(i in colnames(beta_mat)){
    reg_res[,i]<- beta_mat[,colnames(beta_mat)==i]
  }
  
  for(i in colnames(p_mat)){
    reg_res[,i]<- p_mat[,colnames(p_mat)==i]
  }
  
  write.table(reg_res,file = "4_DE_downstream/WGCNA/regression_modules_results.txt", sep = ";", quote = F, row.names = T,col.names = T)
  
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 5, 5));
  
  # Display the correlation values within a heatmap plot
  names(datTraits)[6] <- "Cocaine at death"
  
  pdf("/path/to/RNA/DE_downstream/WGCNA/Heatmap_modules_Expr_filtered.pdf", width = 9, height = 19)
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
  
  # Subset for associated modules only that are i) significantly correlated with CUD and ii) are significantly associated with CUD in the regression analysis
  
  pdf("/path/to/RNA/DE_downstream/WGCNA/Heatmap_modules_Expr_filtered_subset.pdf", width = 9, height = 4)
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
  
  # Define variable weight containing the weight column of datTrait
  CUD = as.data.frame(datTraits$CUD)
  names(CUD) = "CUD"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, CUD, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(CUD), sep="");
  names(GSPvalue) = paste("p.GS.", names(CUD), sep="")
  
  mod <- DF(moduleTraitCor)
  
  module <- rownames(moduleTraitPvalue[unname(moduleTraitPvalue[,1]) < 0.05,])
  module <- substr(module, 3, nchar(module))
  
  for(i in module){
    column = match(i, modNames);
    moduleGenes = moduleColors==i;
    
    png(paste0("/path/to/RNA/DE_downstream/WGCNA/module_membership_gene_significance_prot_",i,".png",sep=""),height=12,width=12,res=300,units="in")
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", i, "module"),
                       ylab = "Gene significance for CUD status",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 3, cex.lab = 3, cex.axis = 3, col = i)
    dev.off()
    
  }
  
write.table(mod, "/path/to/RNA/DE_downstream/WGCNA/mod_trait_cor_Expr.txt", sep = ";", quote = F)  

#### Protein data ####
dir0 <- "/path/to/protein/data_analysis_results_V1/"

setwd(dir0)

#read pheno
pheno <- read.csv("/path/to/protein/data_analysis_results_V1/pheno_merged_with_proteomics.txt", sep=";")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#### Protein ####

# Load the expression and trait data saved in the first part
setwd("WGCNA/")
lnames = load(file ="VS_Prot_dataInput_processed.RData") 
#The variable lnames contains the names of loaded variables.
lnames

# Load network data saved in the second part.
lnames = load(file = "/path/to/protein/VS_networkConstruction_Prot.RData");
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
write.table(MEs, file = "MEs_Prot_filtered.txt", sep = ";", quote = F, row.names = T)

# Regression analysis for the modules showing significant correlation with CUD status
ME_reg_dat <- merge(MEs,datTraits, by=0)
ME_reg_dat$CUD <- as.factor(ME_reg_dat$CUD)
mods <- colnames(MEs)
beta_mat <- matrix(0, nrow=5,ncol=length(mods))
colnames(beta_mat)<-mods
rownames(beta_mat) <- colnames(datTraits)[1:5]
p_mat <- matrix(0, nrow=5,ncol=length(mods))
colnames(p_mat)<-mods
rownames(p_mat) <- colnames(datTraits)[1:5]

for(i in mods){
  lm <- lm(get(i)~CUD+Age+pH+PMI+batch, data = ME_reg_dat)
  beta_mat[,i] <- summary(lm)$coefficients[-1,1]
  p_mat[,i] <-summary(lm)$coefficients[-1,4]
}

p_mat[1,][p_mat[1,] <0.05 ] # modules significantly associated with CUD in the regression analysis
#MEgreenyellow      MEyellow       MEbrown   MElightcyan     MEdarkred 
#0.030650682   0.015476448   0.009717392   0.034602037   0.046830022 

colnames(p_mat)<-paste0("p_",colnames(p_mat))
reg_res <- as.data.frame(matrix(0,nrow=5,ncol=46))
colnames(reg_res)<- c(rbind(colnames(beta_mat),colnames(p_mat)))
rownames(reg_res)<-rownames(beta_mat)


for(i in colnames(beta_mat)){
  reg_res[,i]<- beta_mat[,colnames(beta_mat)==i]
}

for(i in colnames(p_mat)){
  reg_res[,i]<- p_mat[,colnames(p_mat)==i]
}

write.table(reg_res,file = "regression_modules_results.txt", sep = ";", quote = F, row.names = T,col.names = T)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
names(datTraits)[6] <- "Cocaine at death"

pdf("/path/to/protein/data_analysis_results_V1/WGCNA/Heatmap_modules_Protein.pdf", width = 5, height = 8)
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

# Subset for modules only that are i) significantly correlated with CUD and ii) are significantly associated with CUD in the regression analysis

pdf("/path/to/protein/data_analysis_results_V1/WGCNA/Heatmap_modules_Prot_filtered_subset.pdf", width = 6.25, height = 3.5)
par(mar = c(6, 8.5, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor[rownames(moduleTraitCor) %in% c("MEyellow","MEbrown"),c(1,6,2:5)],
               xLabels = names(datTraits)[c(1,6,2:5)],
               yLabels = c("MEyellow","MEbrown"),
               ySymbols = c("MEyellow","MEbrown"),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[c(6,7),c(1,6,2:5)],
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Define variable weight containing the weight column of datTrait
CUD = as.data.frame(datTraits$CUD)
names(CUD) = "CUD"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, CUD, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(CUD), sep="");
names(GSPvalue) = paste("p.GS.", names(CUD), sep="")

mod <- DF(moduleTraitCor)

module <- rownames(moduleTraitPvalue[unname(moduleTraitPvalue[,1]) < 0.05,])
module <- substr(module, 3, nchar(module))

for(i in module){
  column = match(i, modNames);
  moduleGenes = moduleColors==i;
  
  png(paste0("/path/to/protein/data_analysis_results_V1/WGCNA/module_membership_gene_significance_prot_",i,".png",sep=""),height=12,width=12,res=300,units="in")
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", i, "module"),
                     ylab = "Gene significance for CUD status",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 3, cex.lab = 3, cex.axis = 3, col = i)
  dev.off()
  
}

write.table(mod, "/path/to/protein/data_analysis_results_V1/WGCNA/mod_trait_cor_Protein.txt", sep = ";", quote = F)

