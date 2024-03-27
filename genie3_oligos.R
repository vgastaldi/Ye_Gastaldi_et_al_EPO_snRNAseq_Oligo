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
registerDoParallel(11)
nCores = 11
plan("multisession", workers = 11)
options(future.globals.maxSize = 64 * 1024^3,future.seed = T)

#### Loading base RDS file ####
readRDS("paper_lineages5.rds") -> paper_lineages5

#### Creating new column ####
paste(paper_lineages5@meta.data$treatment,"_",paper_lineages5@meta.data$clusters_updated,sep = "") -> paper_lineages5@meta.data$lineages_treatment

#### Starting ####
raw_counts_filt <- as.matrix(GetAssayData(paper_lineages5, slot = "counts", assay = "RNA"))

cellInfo <- as.data.frame(paper_lineages5@meta.data$lineages_treatment)
row.names(cellInfo) <- row.names(paper_lineages5@meta.data)
colnames(cellInfo) <- "lineages_treatment"

cellInfo$seuratCluster <- paste(cellInfo$lineages_treatment,row.names(cellInfo),sep="/")
raw_counts_filt <- raw_counts_filt[rowSums(raw_counts_filt) > 0, ]
raw_counts_filt <- raw_counts_filt[, colSums(raw_counts_filt > 0) >= 1000]
dim(raw_counts_filt)
write.csv(raw_counts_filt,"counts_oligos_filter1000genes.csv",row.names = T,col.names = T,quote = F)

# Split genes
genesSplit <- split(sort(row.names(raw_counts_filt)), 1:length(row.names(raw_counts_filt)))
length(genesSplit)
saveRDS(genesSplit,"geneSplit.rds")

# Parallel running
files <- list.files(pattern = "\\.RData$",path = "GENIE3/")

as.numeric(gsub("GENIE3_weightMatrix_|.RData","",files)) -> files
numbers <- 1:length(genesSplit) #23370
numbers <- numbers[!numbers %in% files]

# Define the number of instances to run
n_instances <- 11

# Define a function that takes a subset of numbers as input and runs the R script
run_script <- function(subset) {
  # Assign the subset of numbers to the 'teste' variable in the global environment
  assign("to_do", subset, envir = .GlobalEnv)
  
  # Run the R script
  source("trial_genie3.R")
}

# Split the numbers vector into subsets
subsets <- split(numbers, ceiling(seq_along(numbers)/(length(numbers)/n_instances)))

# Run the function on each subset in parallel
results <- mclapply(subsets, run_script, mc.cores = n_instances)

# Merge results:
library(GENIE3)
linkList_list <- list()
for(i in 1:23370){
  load(paste0("GENIE3/GENIE3_weightMatrix_",i,".RData"))
  linkList_list[[i]] <- getLinkList(weightMatrix)
}
length(linkList_list)
sapply(linkList_list, nrow)
linkList_list <- Filter(function(x) nrow(x) > 0, linkList_list)

library(dplyr)
linkList_list <- lapply(linkList_list, function(df) {
  df$targetGene <- as.character(df$targetGene)
  df
})

linkList <- bind_rows(linkList_list)

colnames(linkList) <- c("TF", "Target", "weight")
linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
#linkList <- linkList[which(linkList[,"weight"]>0.001),]
nrow(linkList)
head(linkList)
save(linkList, file="GENIE3_linkList.RData")

#### Restart ####
load("GENIE3_linkList.RData")

quantile(linkList$weight, probs=c(0.25, 0.75, 0.90)) # Just checking the distribution
dim(linkList)



head(linkList) # this should look like the given below

#P1	P2	weight
#RPS4Y1 	EIF1AY 	0.07585009
#DDX3Y 	EIF1AY 	0.07349344
#ARMCX2 	OC90 	0.07199847
#KDM5D 	EIF1AY 	0.07177748
#MT1G 	MT1M 	0.06332241
#HSPA1A 	JUN 	0.06312445

# checking the weight of pairwise comparison
#plot(linkList$weight, type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
#     ylab="Weight", xlab="Links sorted decreasingly")
#abline(h=0.001, col="blue") # Threshold
#dev.off()

sum(linkList$weight>0.001)/nrow(linkList) ## checking that what fraction of pairs are positive

linkList_001 <- linkList[which(linkList[,"weight"]>0.001),]
#linkList_001 <- linkList ## assigning to other data so that original files doesn't get messed up


#colnames(linkList_001)[1:2] <- c("TF", "Target") # renaming it

#Create the gene-sets & save:

tfModules <- list()

linkList_001$TF <- as.character(linkList_001$TF)

linkList_001$Target <- as.character(linkList_001$Target)

head(linkList_001)
#### Create TF-modules: 
# 1: Weight > 0.001 (filtered in previous step)
tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))

# 2: Weight > 0.005
llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))
head(llminW)
head(tfModules)
# 3: Top 50 targets for each TF
# ("w001" should be ordered decreasingly by weight)
tfModules[["top50"]] <- lapply(tfModules[["w005"]], function(x) x[1:(min(length(x), 50))])
tfModules$top50 <- tfModules$top50[names(which(lengths(tfModules$top50) >= 20))]


head(tfModules)


tfModules_top50targets <- (tfModules$top50)

capture.output(print(tfModules_top50targets[which(lapply(tfModules_top50targets, length) > 20)]), file = "tf_modules_top50Targets_minimum20Genes.txt")

## Now calculating the AUC significance

library(AUCell)
library(doParallel)
library(doRNG)
library(reshape2)
library(dplyr)
set.seed(0)

cells_rankings <- AUCell_buildRankings(raw_counts_filt)


tfModules_top50targetsFiltered <- tfModules_top50targets[which(lapply(tfModules_top50targets, length) >= 20)]

cells_AUC <- AUCell_calcAUC(tfModules_top50targetsFiltered, cells_rankings)

par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE) 

# that's it
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellInfo), cellInfo[,"lineages_treatment"]) 
cells_AUC <- cells_AUC[onlyNonDuplicatedExtended(rownames(cells_AUC)),]
auc_matrix <- as.matrix(cells_AUC@assays@data$AUC)

# Initialize an empty matrix to store the results
regulonActivity_byCellType <- matrix(nrow = nrow(auc_matrix), ncol = length(cellsPerCluster))
colnames(regulonActivity_byCellType) <- names(cellsPerCluster)
rownames(regulonActivity_byCellType) <- rownames(auc_matrix)

# Calculate the average for each group
for (i in seq_along(cellsPerCluster)) {
  group <- cellsPerCluster[[i]]
  group_matrix <- auc_matrix[, colnames(auc_matrix) %in% group]
  regulonActivity_byCellType[, i] <- rowMeans(group_matrix, na.rm = TRUE)
}

# Print the result
print(regulonActivity_byCellType)

# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

rss <- calcRSS(AUC=getAUC(cells_AUC), cellAnnotation=cellInfo[colnames(cells_AUC), "lineages_treatment"])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

#### Manu's binarization ####
as.data.frame(regulonActivity_byCellType) -> regulon_binary
# Assuming 'df' is your dataframe
row_means <- rowMeans(regulon_binary)
regulon_binary <- sweep(regulon_binary, 1, row_means, "/")
regulon_binary[regulon_binary < 0.7] <- 0
regulon_binary[regulon_binary > 0] <- 1

regulon_binary[,c(10,9,6,7,8,5,4,1,2,3)] -> regulon_binary2

# Load the library
library(ComplexHeatmap)
library(circlize)

# Assuming regulon_binary2 is your dataframe
df <- regulon_binary2

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1:5])) == 1 & length(unique(row[6:10])) == 1 & unique(row[1:5])[1] != unique(row[6:10])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Separate between the states: developing and more mature
#### OPCs e COPs ####
df <- regulon_binary2[,c(1,2,6,7)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1:2])) == 1 & length(unique(row[3:4])) == 1 & unique(row[1:2])[1] != unique(row[3:4])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Anxa2 regulon!
# Open a connection to the file
fileConn <- file("results_regulon_OPC_COP.txt", open = "w")

for (i in 1:nrow(result)){
  # Get the main name
  main_name <- row.names(result)[i]
  
  # Get the associated values with the name
  associated_values <- tfModules_top50targetsFiltered[[main_name]]
  
  # Write the main name to the file
  write(main_name, fileConn)
  
  # Write an empty line to the file
  write("", fileConn)
  
  # Write the associated values to the file, separated by commas
  write(paste(associated_values, collapse = ", "), fileConn)
}

# Close the connection to the file
close(fileConn)

#### MFOLs MOL1 MOL2 ####
df <- regulon_binary2[,c(3,4,5,8,9,10)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1:3])) == 1 & length(unique(row[4:6])) == 1 & unique(row[1:3])[1] != unique(row[4:6])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Nothing!

#### Individual clusters ####
#### OPC ####
df <- regulon_binary2[,c(1,6)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1])) == 1 & length(unique(row[2])) == 1 & unique(row[1])[1] != unique(row[2])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Open a connection to the file
fileConn <- file("results_regulon_OPC.txt", open = "w")

for (i in 1:nrow(result)){
  # Get the main name
  main_name <- row.names(result)[i]
  
  # Get the associated values with the name
  associated_values <- tfModules_top50targetsFiltered[[main_name]]
  
  # Write the main name to the file
  write(main_name, fileConn)
  
  # Write an empty line to the file
  write("", fileConn)
  
  # Write the associated values to the file, separated by commas
  write(paste(associated_values, collapse = ", "), fileConn)
}

# Close the connection to the file
close(fileConn)

#### COP ####
df <- regulon_binary2[,c(2,7)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1])) == 1 & length(unique(row[2])) == 1 & unique(row[1])[1] != unique(row[2])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Open a connection to the file
fileConn <- file("results_regulon_COP.txt", open = "w")

for (i in 1:nrow(result)){
  # Get the main name
  main_name <- row.names(result)[i]
  
  # Get the associated values with the name
  associated_values <- tfModules_top50targetsFiltered[[main_name]]
  
  # Write the main name to the file
  write(main_name, fileConn)
  
  # Write an empty line to the file
  write("", fileConn)
  
  # Write the associated values to the file, separated by commas
  write(paste(associated_values, collapse = ", "), fileConn)
}

# Close the connection to the file
close(fileConn)

#### MFOL ####
df <- regulon_binary2[,c(3,8)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1])) == 1 & length(unique(row[2])) == 1 & unique(row[1])[1] != unique(row[2])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Open a connection to the file
fileConn <- file("results_regulon_MFOL.txt", open = "w")

for (i in 1:nrow(result)){
  # Get the main name
  main_name <- row.names(result)[i]
  
  # Get the associated values with the name
  associated_values <- tfModules_top50targetsFiltered[[main_name]]
  
  # Write the main name to the file
  write(main_name, fileConn)
  
  # Write an empty line to the file
  write("", fileConn)
  
  # Write the associated values to the file, separated by commas
  write(paste(associated_values, collapse = ", "), fileConn)
}

# Close the connection to the file
close(fileConn)

#### MOL1 ####
df <- regulon_binary2[,c(4,9)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1])) == 1 & length(unique(row[2])) == 1 & unique(row[1])[1] != unique(row[2])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Open a connection to the file
fileConn <- file("results_regulon_MOL1.txt", open = "w")

for (i in 1:nrow(result)){
  # Get the main name
  main_name <- row.names(result)[i]
  
  # Get the associated values with the name
  associated_values <- tfModules_top50targetsFiltered[[main_name]]
  
  # Write the main name to the file
  write(main_name, fileConn)
  
  # Write an empty line to the file
  write("", fileConn)
  
  # Write the associated values to the file, separated by commas
  write(paste(associated_values, collapse = ", "), fileConn)
}

# Close the connection to the file
close(fileConn)

#### MOL2 ####
df <- regulon_binary2[,c(5,10)]

# Create a color mapping where 0 values are white and positive values are black
col_fun <- colorRamp2(c(0, max(df)), c("white", "black"))

# Create the heatmap
Heatmap(df, 
        name = "regulon_binary2", 
        col = col_fun,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

# Not a single regulon exclusive to EPO or Placebo
matching_rows <- apply(df, 1, function(row) {
  length(unique(row[1])) == 1 & length(unique(row[2])) == 1 & unique(row[1])[1] != unique(row[2])[1]
})

# Subset df to only include the matching rows
result <- df[matching_rows,]

# Open a connection to the file
fileConn <- file("results_regulon_MOL2.txt", open = "w")

for (i in 1:nrow(result)){
  # Get the main name
  main_name <- row.names(result)[i]
  
  # Get the associated values with the name
  associated_values <- tfModules_top50targetsFiltered[[main_name]]
  
  # Write the main name to the file
  write(main_name, fileConn)
  
  # Write an empty line to the file
  write("", fileConn)
  
  # Write the associated values to the file, separated by commas
  write(paste(associated_values, collapse = ", "), fileConn)
}

# Close the connection to the file
close(fileConn)

#### Saving files ####
saveRDS(regulon_binary,"regulon_binary.rds")
saveRDS(regulonActivity_byCellType,"regulonActivity_byCellType.rds")
saveRDS(cellInfo,"cellInfo.rds")
saveRDS(auc_matrix,"auc_matrix.rds")