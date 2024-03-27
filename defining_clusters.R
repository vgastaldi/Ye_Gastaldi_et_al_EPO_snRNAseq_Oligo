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
readRDS("EPOsnRNA_reduced_labels.rds") -> EPOsnRNA # Object contains a version of the pyramidal paper rds where cell clusters were defined separately and then combined

#### First look at oligos ####
EPOsnRNA[["subsets_Rp_percent"]] <- PercentageFeatureSet(EPOsnRNA, pattern = "^Rp-")
## Removing previous clustering
EPOsnRNA@meta.data[,-c(14,15,21,22,29)] -> EPOsnRNA@meta.data
EPOsnRNA@meta.data -> all_oligos
all_oligos[which(all_oligos$subsets_Mt_percent < 0.005),] -> all_oligos
all_oligos[grepl("Oligo|oligo|OPC|MOL|MFOL",all_oligos$combined_celltypes),] -> all_oligos

#### Getting ready for clustering ####
oligos_EPOsnRNA <- subset(EPOsnRNA, combined_celltypes %in% all_oligos$combined_celltypes)
rm(list = setdiff(ls(), "oligos_EPOsnRNA"))
saveRDS(oligos_EPOsnRNA,"oligos_EPOsnRNA_21062023.rds")
oligos_EPOsnRNA <- NormalizeData(oligos_EPOsnRNA) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent")) %>% RunPCA(verbose = FALSE)
oligos_harmony <- RunHarmony(oligos_EPOsnRNA, group.by.vars = "orig.ident")
rm(oligos_EPOsnRNA)
pct <- oligos_harmony[["pca"]]@stdev / sum(oligos_harmony[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

oligos_harmony <- FindNeighbors(oligos_harmony, reduction = "harmony", dims = 1:co2)
oligos_harmony <- FindClusters(object = oligos_harmony, resolution = 0.5, graph.name="RNA_snn")
oligos_harmony <- RunUMAP(oligos_harmony, reduction = "harmony", dims = 1:co2)
rm(co2,pct)
saveRDS(oligos_harmony, file = "oligos_harmony_newR.rds")

#### Comparing clustering resolutions ####
oligos_harmony <- FindClusters(object = oligos_harmony,
                               resolution = seq(0.1,0.475,0.025))
Idents(oligos_harmony) <- "combined_celltypes"
dimplot_harmony <- DimPlot(oligos_harmony, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none")
dimplot_harmony
oligos_harmony@meta.data$combined_celltypes -> oligos_harmony@meta.data$cells_figure
gsub("V_","",oligos_harmony@meta.data$cells_figure) -> oligos_harmony@meta.data$cells_figure

#### Genes ####
genes=c("Olig2", "Cspg4", "Bcas1", "Apc","Pdgfra","Cnp","Olig1","Bmp4","Plp1","Sox6","Sox10","Mog")

classic_markers <- c("Olig2","Olig1","Cnp","Plp1")
lineage_discovery <- c("Pdgfra","Cspg4","Plp1","Bcas1","Bmp4","Sox6","Apc","Mog")

FeaturePlot(oligos_harmony, 
            reduction = "umap", 
            features = lineage_discovery, 
            order = TRUE,
            min.cutoff = 'q10', 
            label = F)
nebPlp1 <- plot_density(oligos_harmony, "Plp1", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebOlig2 <- plot_density(oligos_harmony, "Olig2", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebOlig1 <- plot_density(oligos_harmony, "Olig1", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebCnp <- plot_density(oligos_harmony, "Cnp", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
wrap_plots(nebOlig1,nebOlig2,nebCnp,nebPlp1)

ggsave(plot = last_plot(), "umap_classic_markers_nebulosa.png", width = 30, height = 20, units = "cm")

nebPdgfra <- plot_density(oligos_harmony, "Pdgfra", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebCspg4 <- plot_density(oligos_harmony, "Cspg4", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebBcas1 <- plot_density(oligos_harmony, "Bcas1", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebBmp4 <- plot_density(oligos_harmony, "Bmp4", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebSox6 <- plot_density(oligos_harmony, "Sox6", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebApc <- plot_density(oligos_harmony, "Apc", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebMog <- plot_density(oligos_harmony, "Mog", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebSox10 <- plot_density(oligos_harmony, "Sox10", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)

wrap_plots(nebPdgfra,nebCspg4,nebPlp1,nebBcas1,nebBmp4,nebSox6,nebApc,nebMog,nebSox10)
ggsave(plot = last_plot(), "umap_lineage_markers_nebulosa.png", width = 30, height = 20, units = "cm")

DefaultAssay(oligos_harmony) <- "RNA"
cds2 <- as.cell_data_set(oligos_harmony)
cds2@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(oligos_harmony[["RNA"]])

#### Running the basic set of functions to get the pseudotime plot ####
cds2 <- cluster_cells(cds2, reduction_method = "UMAP",random_seed = 123)

harmony_mono3_cluster <- plot_cells(cds2, color_cells_by = "cluster", show_trajectory_graph = FALSE,group_label_size = 5)
harmony_mono3_partition <- plot_cells(cds2, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size = 5)
Idents(oligos_harmony) <- "cells_figure"
dimplot_harmony <- DimPlot(oligos_harmony, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
wrap_plots(dimplot_harmony,harmony_mono3_partition)
ggsave(plot = last_plot(), "umap_partition_allclusters.png", width = 30, height = 20, units = "cm")


cds2 <- learn_graph(cds2, use_partition = F, verbose = FALSE)
harmony_trajectory_nopart <- plot_cells(cds2,
                                        color_cells_by = "cluster",
                                        label_groups_by_cluster=T,
                                        label_leaves=FALSE,
                                        label_roots =F,
                                        group_label_size = 5,
                                        label_branch_points=FALSE)

cds2 <- order_cells(cds2, root_cells = colnames(cds2[,clusters(cds2) == 3]))

harmony_pseudotime_nopart <- plot_cells(cds2,
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

wrap_plots(dimplot_harmony,harmony_pseudotime_nopart)
ggsave(plot = last_plot(), "umap_nopartition_pseudotime.png", width = 30, height = 20, units = "cm")

wrap_plots(nebPlp1,harmony_pseudotime_nopart,nebMog)
ggsave(plot = last_plot(), "umap_mature_markers_nebulosa.png", width = 40, height = 20, units = "cm")


cds2 <- learn_graph(cds2, use_partition = T, verbose = FALSE)

harmony_trajectory_part <- plot_cells(cds2,
                                      color_cells_by = "cluster",
                                      label_groups_by_cluster=T,
                                      label_leaves=FALSE,
                                      label_roots =F,
                                      group_label_size = 5,
                                      label_branch_points=FALSE)

cds2 <- order_cells(cds2, root_cells = colnames(cds2[,clusters(cds2) == 3]))

harmony_pseudotime_part1 <- plot_cells(cds2,
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
cds2 <- order_cells(cds2, root_cells = colnames(cds2[,clusters(cds2) == 6]))

harmony_pseudotime_part3 <- plot_cells(cds2,
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
wrap_plots(dimplot_harmony,harmony_pseudotime_part1,harmony_pseudotime_part3)

ggsave(plot = last_plot(), "umap_partitions_pseudotime.png", width = 40, height = 20, units = "cm")

#### Analysis without identified "mixed" cell types ####
readRDS("~/Documents/ownCloud/snRNA/EPOsnRNA_reduced_labels.rds") -> EPOsnRNA
EPOsnRNA[["subsets_Rp_percent"]] <- PercentageFeatureSet(EPOsnRNA, pattern = "^Rp-")
EPOsnRNA@meta.data -> all_oligos
all_oligos[grepl("Oligo|oligo|OPC|MOL|MFOL",all_oligos$combined_celltypes),] -> all_oligos
all_oligos[!grepl("Astrocyte|Pyr|microglia|interneurons|astrocyte",all_oligos$combined_celltypes),] -> all_oligos
unique(all_oligos$combined_celltypes)
all_oligos[all_oligos$subsets_Mt_percent < 0.005,] -> all_oligos
oligos_EPOsnRNA <- subset(EPOsnRNA, combined_celltypes %in% all_oligos$combined_celltypes & subsets_Mt_percent < 0.005)
saveRDS(oligos_EPOsnRNA,"oligos_EPOsnRNA_noMixture_newR.rds")
rm(EPOsnRNA,all_oligos)
oligos_EPOsnRNA <- NormalizeData(oligos_EPOsnRNA,verbose = FALSE) %>% FindVariableFeatures(verbose = FALSE) %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent"),verbose = FALSE) %>% RunPCA(verbose = FALSE) #Disabled to control more
oligos2 <- RunHarmony(oligos_EPOsnRNA, group.by.vars = "orig.ident",verbose = FALSE)
pct <- oligos2[["pca"]]@stdev / sum(oligos2[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)
oligos2 <- FindNeighbors(oligos2, reduction = "harmony", dims = 1:co2,verbose = FALSE)
oligos2 <- FindClusters(object = oligos2, resolution = 0.5, graph.name="RNA_snn",verbose = FALSE)
oligos2 <- RunUMAP(oligos2, reduction = "harmony", dims = 1:co2,verbose = FALSE)
saveRDS(oligos2,"NEW_oligos_harmony_NoMix_newR.rds")
rm(co2,pct,oligos_EPOsnRNA)

#
oligos_harmony_nomix <- oligos2
rm(oligos2)

oligos_harmony_nomix@meta.data$combined_celltypes -> oligos_harmony_nomix@meta.data$cells_figure
gsub("V_","",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
Idents(oligos_harmony_nomix) <- "cells_figure"
DefaultAssay(oligos_harmony_nomix) <- "RNA"

dimplot_nomix <- DimPlot(oligos_harmony_nomix, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
dimplot_nomix

nebulosa_nomix_allgenes <- plot_density(oligos_harmony_nomix, genes, reduction="umap", size=0.6)
ggsave(plot = nebulosa_nomix_allgenes, "nebulosa_nomix_allgenes.png", width = 40, height = 20, units = "cm")

unique(Idents(oligos_harmony_nomix))
oligos_harmony_nomix@meta.data[grepl("Mature",oligos_harmony_nomix$cells_figure),"cells_figure"] <- "MOL"
gsub("Oligodendrocytes_2","Pre-mOLG",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
oligos_harmony_nomix@meta.data[grepl("committed",oligos_harmony_nomix$cells_figure),"cells_figure"] <- "Pre-mOLG"
gsub("committed oligodendrocyte precursors (COP)","Pre-mOLG",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
gsub("MFOL","Possible OL",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
gsub("Oligodendrocytes","Possible OL",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
gsub("OPC_2","Possible OPC",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
gsub("Mature Possible OL (MOL)","MOL",oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
as.character(oligos_harmony_nomix@meta.data$cells_figure) -> oligos_harmony_nomix@meta.data$cells_figure
Idents(oligos_harmony_nomix) <- "cells_figure"
DimPlot(oligos_harmony_nomix, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
ggsave(plot = last_plot(), "umap_oligos_05_nomix.png", width = 30, height = 20, units = "cm")

nebPlp1 <- plot_density(oligos_harmony_nomix, "Plp1", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebOlig2 <- plot_density(oligos_harmony_nomix, "Olig2", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebOlig1 <- plot_density(oligos_harmony_nomix, "Olig1", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebCnp <- plot_density(oligos_harmony_nomix, "Cnp", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
wrap_plots(nebOlig1,nebOlig2,nebCnp,nebPlp1)

ggsave(plot = last_plot(), "umap_classic_markers_nebulosa_nomix.png", width = 30, height = 20, units = "cm")

nebPdgfra <- plot_density(oligos_harmony_nomix, "Pdgfra", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebCspg4 <- plot_density(oligos_harmony_nomix, "Cspg4", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebBcas1 <- plot_density(oligos_harmony_nomix, "Bcas1", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebBmp4 <- plot_density(oligos_harmony_nomix, "Bmp4", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebSox6 <- plot_density(oligos_harmony_nomix, "Sox6", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebApc <- plot_density(oligos_harmony_nomix, "Apc", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebMog <- plot_density(oligos_harmony_nomix, "Mog", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)
nebSox10 <- plot_density(oligos_harmony_nomix, "Sox10", reduction="umap", size=0.6) + scale_color_viridis_c(option="inferno", direction=-1)

wrap_plots(nebPdgfra,nebCspg4,nebPlp1,nebBcas1,nebBmp4,nebSox6,nebApc,nebMog,nebSox10)
ggsave(plot = last_plot(), "umap_lineage_markers_nebulosa_nomix.png", width = 30, height = 20, units = "cm")

#### Monocle3 for reduced cell selection ####
cds_nomix <- as.cell_data_set(oligos_harmony_nomix)
cds_nomix@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(oligos_harmony_nomix[["RNA"]])

cds_nomix <- cluster_cells(cds_nomix, reduction_method = "UMAP",random_seed = 123)

harmony_nomix_mono3_cluster <- plot_cells(cds_nomix, color_cells_by = "cluster", show_trajectory_graph = FALSE,group_label_size = 5)
harmony_nomix_mono3_partition <- plot_cells(cds_nomix, color_cells_by = "partition", show_trajectory_graph = FALSE,group_label_size = 5)
Idents(oligos_harmony_nomix) <- "cells_figure"
dimplot_harmony_nomix <- DimPlot(oligos_harmony_nomix, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
wrap_plots(dimplot_harmony_nomix,harmony_nomix_mono3_cluster,harmony_nomix_mono3_partition)
ggsave(plot = last_plot(), "umap_partition_allclusters_nomix.png", width = 40, height = 20, units = "cm")


cds_nomix <- learn_graph(cds_nomix, use_partition = F, verbose = FALSE)
harmony_nomix_trajectory_nopart <- plot_cells(cds_nomix,
                                              color_cells_by = "cluster",
                                              label_groups_by_cluster=T,
                                              label_leaves=FALSE,
                                              label_roots =F,
                                              group_label_size = 5,
                                              label_branch_points=FALSE)

cds_nomix <- order_cells(cds_nomix, root_cells = colnames(cds_nomix[,clusters(cds_nomix) == 3]))

harmony_nomix_pseudotime_nopart <- plot_cells(cds_nomix,
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

wrap_plots(dimplot_harmony_nomix,harmony_nomix_pseudotime_nopart)
ggsave(plot = last_plot(), "umap_nopartition_pseudotime_nomix.png", width = 30, height = 20, units = "cm")

wrap_plots(nebPlp1,harmony_nomix_pseudotime_nopart,nebMog)
ggsave("umap_mature_markers_nebulosa_nomix.png", width = 40, height = 20, units = "cm")


cds_nomix <- learn_graph(cds_nomix, use_partition = T, verbose = FALSE)

harmony_nomix_trajectory_part <- plot_cells(cds_nomix,
                                            color_cells_by = "cluster",
                                            label_groups_by_cluster=T,
                                            label_leaves=FALSE,
                                            label_roots =F,
                                            group_label_size = 5,
                                            label_branch_points=FALSE)

cds_nomix <- order_cells(cds_nomix, root_cells = colnames(cds_nomix[,clusters(cds_nomix) == 3]))

harmony_nomix_pseudotime_part2 <- plot_cells(cds_nomix,
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
cds_nomix <- order_cells(cds_nomix, root_cells = colnames(cds_nomix[,clusters(cds_nomix) == 6]))

harmony_nomix_pseudotime_part1 <- plot_cells(cds_nomix,
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
wrap_plots(dimplot_harmony_nomix,harmony_nomix_pseudotime_part2,harmony_nomix_pseudotime_part1)

ggsave(plot = last_plot(), "umap_partitions_pseudotime_nomix.png", width = 40, height = 20, units = "cm")

#### The population that branches off has to be explored in another moment ####
## This includes the clusters that were separated in the nomix object
## Subsetting the object
rm(oligos_harmony_nomix)
colnames(oligos_harmony@meta.data)
Idents(oligos_harmony) <- "RNA_snn_res.0.35"
DimPlot(oligos_harmony, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
paper_lineages <- subset(oligos_harmony, RNA_snn_res.0.35 != 11 & RNA_snn_res.0.35 != 10 & RNA_snn_res.0.35 != 9 &
                           RNA_snn_res.0.35 != 8 & RNA_snn_res.0.35 != 7 & RNA_snn_res.0.35 != 5 & RNA_snn_res.0.35 != 4)
DimPlot(paper_lineages, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
plot_density(paper_lineages, genes, reduction="umap", size=0.6)

paper_lineages <- NormalizeData(paper_lineages,verbose = FALSE) %>% FindVariableFeatures(verbose = FALSE) %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent"),verbose = FALSE) %>% RunPCA(verbose = FALSE) #Disabled to control more
paper_lineages2 <- RunHarmony(paper_lineages, group.by.vars = "orig.ident",verbose = FALSE)

pct <- paper_lineages2[["pca"]]@stdev / sum(paper_lineages2[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

paper_lineages2 <- FindNeighbors(paper_lineages2, reduction = "harmony", dims = 1:co2,verbose = FALSE)
paper_lineages2 <- FindClusters(object = paper_lineages2, resolution = 0.5, graph.name="RNA_snn",verbose = FALSE)
paper_lineages2 <- RunUMAP(paper_lineages2, reduction = "harmony", dims = 1:co2,verbose = FALSE)
rm(co2,pct,paper_lineages)
DimPlot(paper_lineages2, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
plot_density(paper_lineages2, genes, reduction="umap", size=0.6)
paper_lineages2 <- subset(paper_lineages2, RNA_snn_res.0.5 != 11)

paper_lineages2 <- NormalizeData(paper_lineages2,verbose = FALSE) %>% FindVariableFeatures(verbose = FALSE) %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent"),verbose = FALSE) %>% RunPCA(verbose = FALSE) #Disabled to control more
paper_lineages3 <- RunHarmony(paper_lineages2, group.by.vars = "orig.ident",verbose = FALSE)

pct <- paper_lineages3[["pca"]]@stdev / sum(paper_lineages3[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

paper_lineages3 <- FindNeighbors(paper_lineages3, reduction = "harmony", dims = 1:co2,verbose = FALSE)
paper_lineages3 <- FindClusters(object = paper_lineages3, resolution = 0.5, graph.name="RNA_snn",verbose = FALSE)
paper_lineages3 <- RunUMAP(paper_lineages3, reduction = "harmony", dims = 1:co2,verbose = FALSE)
rm(co2,pct,paper_lineages2)
DimPlot(paper_lineages3, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
plot_density(paper_lineages3, classic_markers, reduction="umap", size=0.6)
# Some are not lighting up in the MOL clusters, have to increase resolution
paper_lineages3 <- FindClusters(object = paper_lineages3, resolution = seq(0.5,2.5,0.5), graph.name="RNA_snn",verbose = FALSE)
Idents(paper_lineages3) <- "RNA_snn_res.2.5"
DimPlot(paper_lineages3, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")


paper_lineages3 <- subset(paper_lineages3, RNA_snn_res.2.5 != 32)
paper_lineages3 <- NormalizeData(paper_lineages3,verbose = FALSE) %>% FindVariableFeatures(verbose = FALSE) %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent"),verbose = FALSE) %>% RunPCA(verbose = FALSE) #Disabled to control more
paper_lineages4 <- RunHarmony(paper_lineages3, group.by.vars = "orig.ident",verbose = FALSE)

pct <- paper_lineages4[["pca"]]@stdev / sum(paper_lineages4[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

paper_lineages4 <- FindNeighbors(paper_lineages4, reduction = "harmony", dims = 1:co2,verbose = FALSE)
paper_lineages4 <- FindClusters(object = paper_lineages4, resolution = 0.5, graph.name="RNA_snn",verbose = FALSE)
paper_lineages4 <- RunUMAP(paper_lineages4, reduction = "harmony", dims = 1:co2,verbose = FALSE)
rm(co2,pct,paper_lineages3)

DimPlot(paper_lineages4, reduction = "umap",label = T,label.size = 5,repel = T,raster = F) + theme(legend.position="none")
plot_density(paper_lineages4, genes, reduction="umap", size=0.6)

#### Integration with Hou et al dataset ####
## https://doi.org/10.1016/j.celrep.2023.112293
## Trying to further identify clusters
## Nicely provided by the authors
readRDS("/Users/gastaldi/Documents/ownCloud/snRNA/oligos_june/PaperCellReportsColonnaoligo_final v3.rds") -> houV3
#readRDS("~/Downloads/PaperCellReportsColonnaoligo_final v1.rds") -> houV1

gsub("A","C",houV3@meta.data$orig.ident) -> houV3@meta.data$orig.ident
gsub("B","D",houV3@meta.data$orig.ident) -> houV3@meta.data$orig.ident
DefaultAssay(houV3) <- "RNA"
# the name is just to make my life easier, but the definition comes from the paper
as.character(houV3$integrated_snn_res.0.3) -> houV3$combined_SeuratMonocle
unique(houV3$combined_SeuratMonocle)
gsub("1","P_OPC",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("0","P_MOL1",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("2","P_MOL2",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("3","P_DOL",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("4","P_MOL3",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("5","P_MFOL",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("6","P_NFOL/COPs",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle
gsub("7","P_POPC",houV3$combined_SeuratMonocle) -> houV3$combined_SeuratMonocle

# Access the row names of the Seurat object's RNA assay
row_names <- row.names(houV3)

# Use gsub to modify the row names
new_row_names <- gsub("mm10---","",row_names)

# Update the row names of all relevant slots in the Seurat object
rownames(houV3@assays$RNA@counts) <- new_row_names
rownames(houV3@assays$RNA@data) <- new_row_names
rownames(houV3@assays$RNA@scale.data) <- new_row_names

# Combine the two Seurat objects into a single object
combined <- merge(houV3, y = paper_lineages4)
combined[["subsets_Rp_percent"]] <- PercentageFeatureSet(combined, pattern = "^Rp-")
combined[["subsets_Mt_percent"]] <- PercentageFeatureSet(combined, pattern = "^Mt-")

## Adding cell cycle score ##
cell_cycle_genes <- read.delim("/Users/gastaldi/Documents/ownCloud/snRNA/Mus_musculus_cell_cycle_genes.csv",sep=",",row.names = 1)
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Calculate CCS
combined <- CellCycleScoring(combined, g2m.features=g2m_genes, s.features=s_genes)

# Integrate
combined <- NormalizeData(combined) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent")) %>% RunPCA(verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = "orig.ident")

pct <- combined[["pca"]]@stdev / sum(combined[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:co2)
combined <- FindClusters(object = combined, resolution = 0.5, graph.name="RNA_snn")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:co2)
rm(co2,pct)

Idents(combined) <- "combined_SeuratMonocle"
combined_dimplot <- DimPlot(combined, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none")

VlnPlot(combined,c("Pdgfra","Bmp4","Mog","Plp1","Fyn","Olig1","Olig2"))

Idents(combined) <- "combined_SeuratMonocle"
DimPlot(combined, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none")

## Integration worked!
Idents(combined) <- "combined_SeuratMonocle"
## Cleaning all extra clusterization runs
combined@meta.data[,-c(7,9,37:51,53:56)] -> combined@meta.data

## Looking at the one freshly done
Idents(combined) <- "RNA_snn_res.0.5"
DimPlot(combined, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none")
## OPCs are super curious as there is a specific population, cluster 9, that shows up with high markers for Mbp and Plp1
## Cluster 8 is perfect for NFOL/COPs
## Cluster 7 are the MFOL
## Cluster 0 and 1 are the Early MOLs
# Cluster 4 are the Late MOLs
# Cluster 6 is mainly comprised by MOLs coming from the other paper

as.character(combined@meta.data$RNA_snn_res.0.5) -> combined@meta.data$new_labels
gsub("^0","MOL1",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^1","MOL1",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^2","MOL2",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^3","OPC",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^4","MOL2",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^5","OPC",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^6","MOL3",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^7","MFOL",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^8","NFOL/COPs",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^9","OPC2",combined@meta.data$new_labels) -> combined@meta.data$new_labels
gsub("^10","POPC",combined@meta.data$new_labels) -> combined@meta.data$new_labels

# Transfer the labels
paper_lineages4 <- bckp
combined@meta.data[!is.na(combined@meta.data$simple_identities),"new_labels"] -> clusters_updated
names(clusters_updated) <- row.names(combined@meta.data[!is.na(combined@meta.data$simple_identities),])
paper_lineages4[["clusters_updated"]] <- clusters_updated

Idents(paper_lineages4) <- "clusters_updated"
DimPlot(paper_lineages4, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none")
table(paper_lineages4@meta.data$clusters_updated)

VlnPlot(paper_lineages4,c("Pdgfra","Mbp","Cspg4"))
## There is this very weird OPC population, OPC2 with very high Mbp
## It will be removed
## But also OPCs without Pdgfra?
# Identify cells that do not express the gene of interest
cells_noPdgfra <- WhichCells(paper_lineages4, expression = Pdgfra == 0) # use this to calculate number of cells expressing gene

# Create a new Seurat object containing only the selected cells
noPdgfra_pl4 <- subset(paper_lineages4, cells = cells_noPdgfra)
VlnPlot(noPdgfra_pl4,c("Pdgfra","Mbp","Cspg4"))
table(noPdgfra_pl4$clusters_updated)
# Quite a lot of these cells have Mbp

# How many don't have either Pdgfra or Cpsg4?
# First quick of the OPC2 cluster
paper_lineages5 <- subset(paper_lineages4, clusters_updated != "OPC2")

combined <- NormalizeData(paper_lineages5) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress = c("subsets_Mt_percent","G2M.Score","S.Score","subsets_Rp_percent")) %>% RunPCA(verbose = FALSE)
paper_lineages5 <- RunHarmony(paper_lineages5, group.by.vars = "orig.ident")

pct <- paper_lineages5[["pca"]]@stdev / sum(paper_lineages5[["pca"]]@stdev) * 100
co2 <- as.numeric(sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

paper_lineages5 <- FindNeighbors(paper_lineages5, reduction = "harmony", dims = 1:co2)
paper_lineages5 <- FindClusters(object = paper_lineages5, resolution = seq(0.1,0.5,0.05), graph.name="RNA_snn")
paper_lineages5 <- RunUMAP(paper_lineages5, reduction = "harmony", dims = 1:co2)
rm(co2,pct)
Idents(paper_lineages5) <- "clusters_updated"
Idents(paper_lineages5) <- "RNA_snn_res.0.25"
table(paper_lineages5$clusters_updated)
wrap_plots(DimPlot(paper_lineages5, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none"),VlnPlot(paper_lineages5,c("Pdgfra","Mbp","Cspg4")))


# Merge MOL1 and MOL3 as they are mainly located together
DimPlot(paper_lineages5, reduction = "umap",label = T,label.size = 3,repel = T,raster = F)
gsub("MOL3","MOL1",paper_lineages5@meta.data$clusters_updated) -> paper_lineages5@meta.data$clusters_updated
wrap_plots(DimPlot(paper_lineages5, reduction = "umap",label = T,label.size = 3,repel = T,raster = F) + theme(legend.position="none"),VlnPlot(paper_lineages5,c("Pdgfra","Mbp","Cspg4")))

DimPlot(paper_lineages5, reduction = "umap",label = T,label.size = 9,repel = T,raster = F) + theme(legend.position="none")
Idents(paper_lineages5) <- "clusters_plotting"
DimPlot(paper_lineages5, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, raster = FALSE, cols = c("OPC" = "#b38d84","COPs" = "#5d9166", "MFOL" = "#db6060", "MOL1" = "#ed7cbe", "MOL2" = "#5b72bd"),pt.size = 0.3,label.box = T) + theme(legend.position="none")

ggsave(plot = last_plot(), "PAPER_umap_final_idents_paperlineages5.png", width = 14, height = 10, units = "cm",dpi = 600)
saveRDS(paper_lineages5,"paper_lineages5.rds")