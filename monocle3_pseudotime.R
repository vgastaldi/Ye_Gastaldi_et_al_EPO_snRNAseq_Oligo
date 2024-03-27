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

#### monocle 3 ####
cds_pl5 <- as.cell_data_set(paper_lineages5)
cds_pl5@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(paper_lineages5[["RNA"]])

cds_pl5 <- cluster_cells(cds_pl5, reduction_method = "UMAP",random_seed = 0)

harmony_nomix_mono3_cluster <- plot_cells(cds_pl5, color_cells_by = "cluster", show_trajectory_graph = FALSE,group_label_size = 5)
harmony_nomix_mono3_partition <- plot_cells(cds_pl5, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size = 5)
wrap_plots(harmony_nomix_mono3_cluster,harmony_nomix_mono3_partition)
ggsave(plot = last_plot(), "umap_partition_allclusters_nomix_pl5.png", width = 40, height = 20, units = "cm")


cds_pl5 <- learn_graph(cds_pl5, use_partition = F, verbose = FALSE)
harmony_cds_pl5_nopart <- plot_cells(cds_pl5,
                                     color_cells_by = "cluster",
                                     label_groups_by_cluster=T,
                                     label_leaves=FALSE,
                                     label_roots =F,
                                     group_label_size = 5,
                                     label_branch_points=FALSE)

cds_pl5 <- order_cells(cds_pl5, root_cells = colnames(cds_pl5[,clusters(cds_pl5) == 3]))

harmony_cds_pl5_pseudotime_nopart <- plot_cells(cds_pl5,
                                                color_cells_by = "pseudotime",
                                                group_cells_by = "cluster",
                                                label_cell_groups = FALSE,
                                                label_groups_by_cluster=T,
                                                label_leaves=F,
                                                label_branch_points=F,
                                                label_roots = F,
                                                trajectory_graph_color = "grey28",trajectory_graph_segment_size = 1) + theme(legend.key.size = unit(0.25,"cm")) + theme(
                                                  legend.title = element_text(size = 5),
                                                  legend.text = element_text(size = 4)
                                                ) + scale_colour_gradient2(low = "lightblue", mid = "gold", high = "blue")
harmony_cds_pl5_pseudotime_nopart
ggsave("harmony_cds_pl5_pseudotime_nopart.png",harmony_cds_pl5_pseudotime_nopart,width = 25,height = 17,units = "cm",dpi = 600)

wrap_plots(DimPlot(paper_lineages5, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none"),harmony_cds_pl5_pseudotime_nopart) -> cluster_and_pseudotime

ggsave("cluster_and_pseudotime_paperlineages5.png",cluster_and_pseudotime,width = 25,height = 17,units = "cm",dpi = 600)

#### Pseudotime comparison ####
## First from monocle3 ##
paper_lineages5[["pseudotime_monocle3"]] <- pseudotime(paperL5_monocle)

## Also getting it from monocle2 DDRTree ##
data <- as(as.matrix(paper_lineages5@assays$RNA@data), 'sparseMatrix')
DefaultAssay(paper_lineages5) = "RNA"

pd <- new('AnnotatedDataFrame', data = paper_lineages5@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 1,
                       expressionFamily = uninormal())
gc()

#### Reducing dimension through DDRTree 
HSMM <- reduceDimension(HSMM,norm_method="none", reduction_method="DDRTree",  max_components=4, scaling=TRUE,    verbose=TRUE,  pseudo_expr=1)
HSMM <- orderCells(HSMM)

#### Comparing pseudotime of treatments ####
HSMM@phenoData@data$Pseudotime -> pseudotime_monocle2
colnames(HSMM) -> names(pseudotime_monocle2)
paper_lineages5[["pseudotime_monocle2"]] <- pseudotime_monocle2 # makes absolute no sense in the plot 
boxplot_pseudotime <- ggplot(paper_lineages5@meta.data,aes(x=clusters_updated,y=pseudotime_monocle3,fill=treatment),color = c("blue","red"))+ geom_boxplot(outlier.shape = NA) + stat_summary(fun=mean, geom="point", shape=20, size=3, color="grey10",position = position_dodge2(width = 0.75,preserve = "single")) + scale_fill_manual(values = c("royalblue2","red2")) + theme_minimal() 
boxplot_pseudotime

pseudotime_per_group$treatment <- NA
pseudotime_per_group[pseudotime_per_group$type == "A_Samples","treatment"] <- "Epo"
pseudotime_per_group[pseudotime_per_group$type == "B_Samples","treatment"] <- "Placebo"
as.character(pseudotime_per_group$celltypes) -> pseudotime_per_group$celltypes_chr
as.factor(pseudotime_per_group$treatment) -> pseudotime_per_group$treatment
#as.factor(pseudotime_per_group$celltypes_chr) -> pseudotime_per_group$celltypes_chr
levels(pseudotime_per_group$treatment) <- c("Placebo","Epo")
pseudotime_per_group$celltypes <- factor(pseudotime_per_group$celltypes,levels=c("New_born_migrating_superficial_Sox5","New_born_migrating_superficial","Ventral_CA1_migrating","New_born_migrating_serotonin_firing",
                                                                                 "CA1_superficial","CA1_dorsal"))
boxplot_pseudotime <- ggplot(pseudotime_per_group,aes(x=celltypes,y=Pseudotime,fill=treatment),color = c("blue","red"))+ geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="grey10",position = position_dodge2(width = 0.75,preserve = "single")) + scale_fill_manual(values = c("royalblue2","red2")) + theme_minimal()
ggsave("boxplot_pseudotime.png",boxplot_pseudotime,width = 25,height = 17,units = "cm",dpi = 600)
# No removal of outliers

wilcox_cluster <- NULL
for (i in unique(paper_lineages5@meta.data$clusters_updated)){
  paper_lineages5@meta.data[paper_lineages5@meta.data$clusters_updated == i,] -> testing_group
  as.data.frame(i) -> wilcox_cluster_hold
  wilcox_cluster_hold$p_value <- wilcox.test(testing_group$pseudotime_monocle3 ~ testing_group$treatment,correct = T)$p.value
  colnames(wilcox_cluster_hold) <- "Cluster"
  rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
  for (x in 1:999){
    sample_n(testing_group, 1000,replace = T) -> teste
    as.data.frame(i) -> wilcox_cluster_hold
    wilcox_cluster_hold$p_value <- wilcox.test(teste$pseudotime_monocle3 ~ teste$treatment,correct = T)$p.value
    colnames(wilcox_cluster_hold) <- "Cluster"
    rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
    x <- x + 1
  }
}
fdr_wilcox_cluster <- NULL
for (i in unique(wilcox_cluster$Cluster)){
  wilcox_cluster[wilcox_cluster$Cluster == i,] -> hold
  p.adjust(hold[,2],"fdr",1000) -> hold
  as.data.frame(i) -> fdr_wilcox_hold
  fdr_wilcox_hold$fdr <- (1000-length(hold[hold <= 0.05]))/1000
  rbind(fdr_wilcox_cluster,fdr_wilcox_hold) -> fdr_wilcox_cluster
}

write.table(fdr_wilcox_cluster,"fdr_wilcox_cluster_pseudotime.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#### EPO vs Placebo - all together ####
wilcox_cluster <- NULL
paper_lineages5@meta.data -> testing_group
as.data.frame("EPO vs Placebo all") -> wilcox_cluster_hold
wilcox_cluster_hold$p_value <- wilcox.test(testing_group$pseudotime_monocle3 ~ testing_group$treatment,correct = T)$p.value
colnames(wilcox_cluster_hold) <- "Cluster"
rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
for (x in 1:999){
  sample_n(testing_group, 1000,replace = T) -> teste
  as.data.frame(i) -> wilcox_cluster_hold
  wilcox_cluster_hold$p_value <- wilcox.test(teste$pseudotime_monocle3 ~ teste$treatment,correct = T)$p.value
  colnames(wilcox_cluster_hold) <- "Cluster"
  rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
  x <- x + 1
}

fdr_wilcox_cluster <- NULL
for (i in unique(wilcox_cluster$Cluster)){
  wilcox_cluster[wilcox_cluster$Cluster == i,] -> hold
  p.adjust(hold[,2],"fdr",1000) -> hold
  as.data.frame(i) -> fdr_wilcox_hold
  fdr_wilcox_hold$fdr <- (1000-length(hold[hold <= 0.05]))/1000
  rbind(fdr_wilcox_cluster,fdr_wilcox_hold) -> fdr_wilcox_cluster
}

#### Comparing between the main groups ####
wilcox_cluster <- NULL
paper_lineages5@meta.data[which(paper_lineages5@meta.data$clusters_updated == "MOL2" | paper_lineages5@meta.data$clusters_updated == "MOL1"),] -> testing_group
as.data.frame("MOL2/MOL1") -> wilcox_cluster_hold
wilcox_cluster_hold$p_value <- wilcox.test(testing_group$pseudotime_monocle3 ~ testing_group$clusters_updated,correct = T)$p.value
colnames(wilcox_cluster_hold) <- "Cluster"
rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
for (x in 1:999){
  sample_n(testing_group, 1000,replace = T) -> teste
  as.data.frame("MOL2/MOL1") -> wilcox_cluster_hold
  wilcox_cluster_hold$p_value <- wilcox.test(teste$pseudotime_monocle3 ~ teste$clusters_updated,correct = T)$p.value
  colnames(wilcox_cluster_hold) <- "Cluster"
  rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
  x <- x + 1
}
paper_lineages5@meta.data[which(paper_lineages5@meta.data$clusters_updated == "MOL1" | paper_lineages5@meta.data$clusters_updated == "MFOL"),] -> testing_group
as.data.frame("MOL1/MFOL") -> wilcox_cluster_hold
wilcox_cluster_hold$p_value <- wilcox.test(testing_group$pseudotime_monocle3 ~ testing_group$clusters_updated,correct = T)$p.value
colnames(wilcox_cluster_hold) <- "Cluster"
rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
for (x in 1:999){
  sample_n(testing_group, 1000,replace = T) -> teste
  as.data.frame("MOL1/MFOL") -> wilcox_cluster_hold
  wilcox_cluster_hold$p_value <- wilcox.test(teste$pseudotime_monocle3 ~ teste$clusters_updated,correct = T)$p.value
  colnames(wilcox_cluster_hold) <- "Cluster"
  rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
  x <- x + 1
}
paper_lineages5@meta.data[which(paper_lineages5@meta.data$clusters_updated == "MFOL" | paper_lineages5@meta.data$clusters_updated == "NFOL/COPs"),] -> testing_group
as.data.frame("MFOL/NFOL_COPs") -> wilcox_cluster_hold
wilcox_cluster_hold$p_value <- wilcox.test(testing_group$pseudotime_monocle3 ~ testing_group$clusters_updated,correct = T)$p.value
colnames(wilcox_cluster_hold) <- "Cluster"
rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
for (x in 1:999){
  sample_n(testing_group, 1000,replace = T) -> teste
  as.data.frame("MFOL/NFOL_COPs") -> wilcox_cluster_hold
  wilcox_cluster_hold$p_value <- wilcox.test(teste$pseudotime_monocle3 ~ teste$clusters_updated,correct = T)$p.value
  colnames(wilcox_cluster_hold) <- "Cluster"
  rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
  x <- x + 1
}
paper_lineages5@meta.data[which(paper_lineages5@meta.data$clusters_updated == "NFOL/COPs" | paper_lineages5@meta.data$clusters_updated == "OPC"),] -> testing_group
as.data.frame("NFOL_COPs/OPC") -> wilcox_cluster_hold
wilcox_cluster_hold$p_value <- wilcox.test(testing_group$pseudotime_monocle3 ~ testing_group$clusters_updated,correct = T)$p.value
colnames(wilcox_cluster_hold) <- "Cluster"
rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
for (x in 1:999){
  sample_n(testing_group, 1000,replace = T) -> teste
  as.data.frame("NFOL_COPs/OPC") -> wilcox_cluster_hold
  wilcox_cluster_hold$p_value <- wilcox.test(teste$pseudotime_monocle3 ~ teste$clusters_updated,correct = T)$p.value
  colnames(wilcox_cluster_hold) <- "Cluster"
  rbind(wilcox_cluster,wilcox_cluster_hold) -> wilcox_cluster
  x <- x + 1
}
fdr_wilcox_cluster <- NULL
for (i in unique(wilcox_cluster$Cluster)){
  wilcox_cluster[wilcox_cluster$Cluster == i,] -> hold
  p.adjust(hold[,2],"fdr",1000) -> hold
  as.data.frame(i) -> fdr_wilcox_hold
  fdr_wilcox_hold$fdr <- (1000-length(hold[hold <= 0.05]))/1000
  rbind(fdr_wilcox_cluster,fdr_wilcox_hold) -> fdr_wilcox_cluster
}
fdr_wilcox_cluster[,2] <- ifelse(fdr_wilcox_cluster[,2] == 0, "< 0.001", fdr_wilcox_cluster[,2])
write.table(fdr_wilcox_cluster,"fdr_wilcox_clusterVScluster_pseudotime.txt",col.names = T,row.names = F,quote = F,sep = "\t")