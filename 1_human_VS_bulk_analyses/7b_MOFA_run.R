## Run MOFA - adapted from Tutorial code https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html and https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html
# Authors: Lea Zillich and Eric Zillich 
# last change: EZ 2024-09-12

library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(psych)
library(ggpubr)

DF <- data.frame

dir0 <- "/path/to/protein/data_analysis_results_V1/MOFA"

setwd(dir0)

#Specify the python version that will be used for training the model
reticulate::use_python("/home/user/anaconda3/bin/python", required=TRUE)

# import input data
VS <- get(load("VS_mRNA_miRNA_protein.Rdata"))

# Add names and bring into the same order
names(VS) <- c("miRNA","mRNA","protein")
VS$miRNA <- VS$miRNA[,order(as.numeric(colnames(VS$miRNA)),decreasing = F)]
VS$mRNA <- VS$mRNA[,order(as.numeric(colnames(VS$mRNA)),decreasing = F)]
VS$protein <- VS$protein[,order(as.numeric(colnames(VS$protein)),decreasing = F)]

  # Create MOFA object
  MOFAobject <- create_mofa(VS)
  
  # Specify MOFA options
  data_opts <- get_default_data_options(MOFAobject)

  model_opts <- get_default_model_options(MOFAobject)
  
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- 42
  train_opts$maxiter <- 10000
  train_opts$drop_factor_threshold <- 0.01

 
  # Train the model
  MOFAobject <- prepare_mofa(MOFAobject,
                             data_options = data_opts,
                             model_options = model_opts,
                             training_options = train_opts
  )
  
  MOFAobject <- run_mofa(MOFAobject, outfile="output/MOFA_VS_trained_model.hdf5",use_basilisk = F)
  saveRDS(MOFAobject,"output/MOFA_VS_trained_model_20240303.rds")



