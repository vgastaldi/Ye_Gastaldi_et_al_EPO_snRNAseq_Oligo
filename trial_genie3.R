#### Oligodendrocytes analysis - SCENIC/GENIE3 ####
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(doParallel)
library(future)
library(cowplot)
library(patchwork)
library(SeuratWrappers)
library(Nebulosa)
library(dplyr)
library(SCENIC)
library(GENIE3)
library(doRNG)
library(SCopeLoomR)
set.seed(0)
registerDoParallel(1)
nCores = 1
plan("multisession", workers = 1)
options(future.globals.maxSize = 64 * 1024^3,future.seed = T)
setwd("~/Documents/ownCloud/snRNA/oligos_june/")

readRDS("geneSplit.rds") -> genesSplit
as.matrix(read.delim("counts_oligos_filter1000genes.csv",sep=",",row.names = 1)) -> raw_counts_filt

# Run code
for(i in to_do){
  print(i)
  set.seed(0)
  weightMatrix <- GENIE3(raw_counts_filt, regulators = NULL,nCores=1, targets=genesSplit[[i]])
  save(weightMatrix, file=paste0("GENIE3/GENIE3_weightMatrix_",i,".RData"))
}