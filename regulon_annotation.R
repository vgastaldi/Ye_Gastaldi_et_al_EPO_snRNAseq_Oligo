# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")  # Mouse genome annotation package
set.seed(0)

# Load packages
library(clusterProfiler)
library(org.Mm.eg.db)

arqvo <- dir(path = getwd(), pattern ="results_regulon")
arqvo[-c(1,6)] -> arqvo

#### Annotation ####

table_GOs <- NULL
table_KEGGs <- NULL

for (i in 1:length(arqvo)){
  # Read the file
  lines <- readLines(arqvo[i])
  # Remove empty lines
  lines <- lines[nzchar(lines)] 
  
  # Initialize an empty list
  data_list <- list()
  
  # Initialize a variable to hold the current name
  current_name <- NULL
  
  # Loop over the lines
  for (line in lines) {
    # Split the line into words
    words <- strsplit(line, ",")[[1]]
    
    # Check if the line contains more than one word
    if (length(words) > 1) {
      # Add the current name to the beginning of the words
      words <- c(current_name, words)
      
      # Add the words to the list under the current name
      data_list[[current_name]] <- words
    } else {
      # Update the current name
      current_name <- words
    }
  }
  for (x in 1:length(data_list)){
    genes <- unlist(data_list[x])
    
    genes[!grepl("^Gm|Rik$|^BC[0-9]",genes)] -> genes
    
    entrez_ids <- mapIds(org.Mm.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL")
    
    ego <- enrichGO(gene         = entrez_ids,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "ALL",  # Biological Process
                    pAdjustMethod = "BH",  # Adjust method
                    minGSSize = 3,
                    pvalueCutoff  = 0.05,  # Adjust p-value cutoff
                    qvalueCutoff  = 0.2,  # Adjust q-value cutoff
                    readable      = TRUE)  # Make the result readable
    
    ego@result -> ego_results
    
    if(nrow(ego_results) > 0){
      ego_results$gene_regulon <- data_list[[x]][1]
      
      ego_results$regulon <- gsub("results_regulon_|.txt","",arqvo[i])
      
      ego_results[,c(ncol(ego_results),(ncol(ego_results)-1),1:(ncol(ego_results)-2))] -> ego_results
      
      rbind(table_GOs,ego_results) -> table_GOs
   
    }
    
    kegg_results <- enrichKEGG(gene         = entrez_ids,
                               organism     = "mmu",  # Mouse
                               keyType = "kegg",
                               pvalueCutoff = 0.05,  # p-value cutoff
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.2, # q-value cutoff
    )
    
    kegg_results_df <- kegg_results@result
    
    kegg_results_df[which(kegg_results_df$p.adjust <= 0.05),] -> kegg_results_df
    
    if(nrow(kegg_results_df) > 0){
      # Convert Entrez IDs to gene symbols
      for (z in 1:nrow(kegg_results_df)) {
        # Get the geneID
        geneID <- kegg_results_df$geneID[z]
        
        # Split the geneID into individual Entrez IDs
        entrez_ids <- strsplit(geneID, "/")[[1]]
        
        # Convert Entrez IDs to gene symbols
        gene_symbols <- select(org.Mm.eg.db, keys=entrez_ids, columns="SYMBOL", keytype="ENTREZID")
        
        # Combine gene symbols into a single string
        gene_symbols_string <- paste(gene_symbols$SYMBOL, collapse="/")
        
        # Replace kegg_results_df$geneID with gene symbols
        kegg_results_df$geneID[z] <- gene_symbols_string
      }
      
      kegg_results_df$gene_regulon <- data_list[[x]][1]
      
      kegg_results_df$regulon <- gsub("results_regulon_|.txt","",arqvo[i])
      
      kegg_results_df[,c(ncol(kegg_results_df),(ncol(kegg_results_df)-1),1:(ncol(kegg_results_df)-2))] -> kegg_results_df
      
      rbind(table_KEGGs,kegg_results_df) -> table_KEGGs
    }
  }
}

openxlsx::write.xlsx(table_GOs,"table_GOs_all_regulons_masterGene.xlsx",rowNames = F)
openxlsx::write.xlsx(table_KEGGs,"table_KEGGs_all_regulons_masterGene.xlsx",rowNames = F)

#### Regulon table ####
regulon_table <- NULL

for (i in 1:length(arqvo)){
  # Read the file
  lines <- readLines(arqvo[i])
  # Remove empty lines
  lines <- lines[nzchar(lines)] 
  
  # Initialize an empty list
  data_list <- list()
  
  # Initialize a variable to hold the current name
  current_name <- NULL
  
  # Loop over the lines
  for (line in lines) {
    # Split the line into words
    words <- strsplit(line, ",")[[1]]
    
    # Check if the line contains more than one word
    if (length(words) > 1) {
      # Add the current name to the beginning of the words
      words <- c(current_name, words)
      
      # Add the words to the list under the current name
      data_list[[current_name]] <- words
    } else {
      # Update the current name
      current_name <- words
    }
  }
  as.data.frame(gsub("results_regulon_|.txt","",arqvo[i])) -> hold
  colnames(hold) <- "cluster"
  for (z in 1:length(data_list)){
    hold -> hold2
    hold2$master_gene <- data_list[[z]][1]
    hold2$genes <- paste(data_list[[z]][-1], collapse = ", ")
    rbind(regulon_table,hold2) -> regulon_table
  }
}

paste(data_list[[z]][-1], collapse = ", ")

library(dplyr)

levels <- c("OPC", "COP", "MFOL", "MOL1", "MOL2")

regulon_table$cluster <- factor(regulon_table$cluster, levels = levels)

# Order the dataframe first by 'regulon' (in the order specified by 'levels') and then by 'gene_regulon' alphabetically
regulon_table <- regulon_table %>%
  arrange(cluster, master_gene)

openxlsx::write.xlsx(regulon_table,"regulon_table.xlsx",rowNames = F)