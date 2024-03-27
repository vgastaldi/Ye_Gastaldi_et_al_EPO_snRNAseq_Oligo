#### Running DGE ####
#### Oligodendrocytes analysis ####
### Pre work ####
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(doParallel)
library(future)
library(cowplot)
library(patchwork)
library(monocle3)
library(monocle)
library(SeuratWrappers)
library(Nebulosa)
library(dplyr)
library(edgeR)
library(harmony)
library(AnnotationHub)
set.seed(0)
registerDoParallel(11)
nCores = 11
plan("multisession", workers = 11)
options(future.globals.maxSize = 56 * 1024^3,future.seed = T)

#### Reading initial file ####
readRDS("paper_lineages5.rds") -> paper_lineages5

DEG_simple <- NULL
for (i in unique(paper_lineages5$clusters_updated)){
  ols <- subset(paper_lineages5, clusters_updated == i & ident != "A5" & ident != "B12")
  Idents(ols) <- "treatment"
  OLSmarkers <- FindMarkers(object = ols, ident.1 = "Erythropoietin", ident.2 = "Placebo",test.use = "bimod",logfc.threshold = 0.15,min.pct = 0.1,random.seed = 0,densify = T)
  OLSmarkers[which(OLSmarkers$p_val_adj < 0.05),] -> OLSmarkers
  row.names(OLSmarkers) -> OLSmarkers$gene
  i -> OLSmarkers$cluster
  rbind(DEG_simple,OLSmarkers) -> DEG_simple
  rm(ols,OLSmarkers)
}
write.table(DEG_simple,"DEG_simple_OligoLineagesEpo.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(DEG_simple[which(DEG_simple$avg_log2FC >= 0.25 | DEG_simple$avg_log2FC <= -0.25),],"DEG_simple_OligoLineagesEpo_FC025.txt",row.names = F,col.names = T,quote = F,sep = "\t")
write.table(DEG_simple[which((DEG_simple$avg_log2FC >= 0.25 | DEG_simple$avg_log2FC <= -0.25) & DEG_simple$cluster == "OPC"),"gene"],"DEG_simple_OligoLineagesEpo_FC025_genesOPC.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(DEG_simple[which((DEG_simple$avg_log2FC >= 0.25 | DEG_simple$avg_log2FC <= -0.25) & DEG_simple$cluster == "NFOL/COPs"),"gene"],"DEG_simple_OligoLineagesEpo_FC025_genesNFOL.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(DEG_simple[which((DEG_simple$avg_log2FC >= 0.25 | DEG_simple$avg_log2FC <= -0.25) & DEG_simple$cluster == "MOL1"),"gene"],"DEG_simple_OligoLineagesEpo_FC025_genesMOL1.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(DEG_simple[which((DEG_simple$avg_log2FC >= 0.25 | DEG_simple$avg_log2FC <= -0.25) & DEG_simple$cluster == "MOL2"),"gene"],"DEG_simple_OligoLineagesEpo_FC025_genesMOL2.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(DEG_simple[which((DEG_simple$avg_log2FC >= 0.25 | DEG_simple$avg_log2FC <= -0.25) & DEG_simple$cluster == "MFOL"),"gene"],"DEG_simple_OligoLineagesEpo_FC025_genesMFOL.txt",row.names = F,col.names = F,quote = F,sep = "\t")