library(Seurat)
library(dplyr)
library(BiocParallel)
library(arrow)
library(SeuratWrappers)
library(ggrepel)
library(future)
library(data.table)
library(tidyverse)
options(future.globals.maxSize = 300 * 1024^3) 
library(dbscan)

################# 1. Perform spatially resolved homotypic immune cell colony analysis
#################
#################
#################
#################
#################

################################################################################
# 1_1. Load data
################################################################################
metadata_merge<-readRDS("~/space/seurat_niche/metadata_merge.rds") # Merged meta data of all cancers was generated in Step02, line 139
table(metadata_merge$celltype)

################################################################################
# 1_2. Select cell types for cell colony analysis
################################################################################
metadata_merge_select<-metadata_merge
celltype_immune<-c("T_cell","B_cell","Plasma","Myeloid","Fibroblast","Smooth_muscle","Endothelial","Pericyte")
metadata_merge_select<-metadata_merge_select[which(metadata_merge_select$celltype %in% celltype_immune),]
dim(metadata_merge_select)
spatial_data<-droplevels(metadata_merge_select)

################################################################################
# 1_3. Define a function to identify cell colony
################################################################################

# Identify colonies for each cell type
identify_colonies <- function(cell_type_data, eps = 30, min_pts = 3, cell_type_name = NULL, cancer_type = NULL) {
  if(nrow(cell_type_data) < min_pts) {
    # If there are too few cells, create independent colony_ids for each cell
    cell_type_data$colony_id <- paste0("single_", seq_len(nrow(cell_type_data)))
    cell_type_data$colony_size <- 1
    
    # Add unique colony identifiers
    if (!is.null(cell_type_name) && !is.null(cancer_type)) {
      cell_type_data$unique_colony_id <- paste0(cancer_type, "_", cell_type_name, "_single", seq_len(nrow(cell_type_data)))
    } else {
      cell_type_data$unique_colony_id <- paste0("single", seq_len(nrow(cell_type_data)))
    }
    
    return(cell_type_data)
  }
  
  # Use DBSCAN clustering to identify colonies
  coords <- cell_type_data[, c("coord_x", "coord_y")]
  dbscan_result <- dbscan(coords, eps = eps, minPts = min_pts)
  
  # Process DBSCAN results, create unique IDs for each single cell
  n_single_cells <- sum(dbscan_result$cluster == 0)
  
  # Create new colony_id
  colony_ids <- character(nrow(cell_type_data))
  
  # Process single cells
  single_mask <- dbscan_result$cluster == 0
  if (n_single_cells > 0) {
    colony_ids[single_mask] <- paste0("single_", seq_len(n_single_cells))
  }
  
  # Process colony cells
  colony_mask <- dbscan_result$cluster > 0
  if (sum(colony_mask) > 0) {
    colony_ids[colony_mask] <- paste0("colony_", dbscan_result$cluster[colony_mask])
  }
  
  cell_type_data$colony_id <- colony_ids
  cell_type_data$colony_size <- ave(rep(1, nrow(cell_type_data)), 
                                    colony_ids, FUN = sum)
  
  # Create unique colony identifiers
  if (!is.null(cell_type_name) && !is.null(cancer_type)) {
    cell_type_data$unique_colony_id <- paste0(cancer_type, "_", cell_type_name, "_", colony_ids)
  } else {
    cell_type_data$unique_colony_id <- colony_ids
  }
  
  return(cell_type_data)
}

# Classify colonies by colony size
classify_colony_size <- function(colony_data) {
  if(nrow(colony_data) == 0) {
    return(data.frame())
  }
  
  # Ensure breaks and labels match
  colony_data$size_group <- cut(colony_data$colony_size,
                                breaks = c(0, 1, 5, 20, 50, Inf),  # Corrected to 5 instead of 10
                                labels = c("single(1)", "small(2-5)", "medium(6-20)", 
                                           "large(21-50)", "huge(>50)"),
                                include.lowest = TRUE)
  
  # For single cells (colony_id starting with "single_"), force set to "single(1)"
  colony_data$size_group[grepl("^single_", colony_data$colony_id)] <- "single(1)"
  
  return(colony_data)
}
################################################################################
# 1_4. Identify cell colonies for each cancer type
################################################################################

cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")

# Create an empty list to store results for each cancer type
all_cancer_results <- list()

for (cancer in cancer_types) {
  print(paste("Processing cancer type:", cancer))
  
  # Select data for the current cancer type
  spatial_data_select <- spatial_data[which(spatial_data$cancertype == cancer), ]
  
  # Check if data exists
  if (nrow(spatial_data_select) == 0) {
    print(paste("No data found for cancer type:", cancer))
    next
  }
  
  colony_results <- list()
  unique_cell_types <- names(table(spatial_data_select$celltype))
  unique_cell_types <- unique_cell_types[!is.na(unique_cell_types)]  # Remove NA
  
  for(cell_type in unique_cell_types) {
    cell_type_data <- spatial_data_select[spatial_data_select$celltype == cell_type, ]
    
    # Check if data exists for this cell type
    if (nrow(cell_type_data) == 0) next
    
    # Pass cell type and cancer type info to generate unique identifiers
    colony_data <- identify_colonies(cell_type_data, cell_type_name = cell_type, cancer_type = cancer)
    colony_results[[cell_type]] <- colony_data
  }
  
  # Remove null results
  colony_results <- colony_results[sapply(colony_results, function(x) !is.null(x) && nrow(x) > 0)]
  
  if (length(colony_results) > 0) {
    # Apply classification for the current cancer type
    classified_colonies <- lapply(colony_results, classify_colony_size)
    
    # Remove empty data frames
    classified_colonies <- classified_colonies[sapply(classified_colonies, nrow) > 0]
    
    if (length(classified_colonies) > 0) {
      # Merge results for current cancer type and add cancer type identifier
      cancer_result <- do.call(rbind, classified_colonies)
      cancer_result$analysis_cancertype <- cancer  # Add label for cancer type used in analysis
      
      # Store in the master result list
      all_cancer_results[[cancer]] <- cancer_result
      print(paste("Processed", cancer, ":", nrow(cancer_result), "cells"))
      
      # Print colony statistics for the current cancer type
      colony_stats <- cancer_result %>%
        group_by(celltype, unique_colony_id, size_group) %>%
        summarise(n_cells = n(), .groups = 'drop')
      
      print(paste("Colony statistics for", cancer, ":"))
      print(table(colony_stats$size_group))
    }
  }
}
table(cancer_result[which(cancer_result$celltype=="T_cell"),]$size_group)

################################################################################
# 1_5. Merge data
################################################################################

if (length(all_cancer_results) > 0) {
  final_result <- do.call(rbind, all_cancer_results)
  
  # Reset row names
  rownames(final_result) <- NULL
  
  print(paste("Final combined result dimensions:", dim(final_result)))
  print("Cancer types processed:")
  print(table(final_result$analysis_cancertype))
  
  # Print final summary of colony information
  print("=== COLONY INFORMATION SUMMARY ===")
  
  # Overall colony 
  total_colonies <- final_result %>%
    filter(colony_id > 0) %>%  # Exclude single cells
    group_by(analysis_cancertype, celltype, unique_colony_id) %>%
    summarise(colony_size = dplyr::first(colony_size),
              size_group = dplyr::first(size_group),
              .groups = 'drop')
  
  print("Total colonies by cancer type and cell type:")
  colony_summary <- total_colonies %>%
    group_by(analysis_cancertype, celltype, size_group) %>%
    summarise(n_colonies = n_distinct(unique_colony_id),
              mean_size = mean(colony_size),
              .groups = 'drop')
  
  print(colony_summary)
  
  # Single cell 
  single_cells <- final_result %>%
    filter(colony_id == 0) %>%
    group_by(analysis_cancertype, celltype) %>%
    summarise(n_single_cells = n(), .groups = 'drop')
  
  print("Single cells by cancer type and cell type:")
  print(single_cells)
} else {
  print("No data processed for any cancer type.")
  final_result <- data.frame()
}
dir.create("~/space/analysis/colony/",recursive = T)
write.csv(final_result,"~/space/analysis/colony/2_4_final_result.csv")
final_result<-read.csv("~/space/analysis/colony/2_4_final_result.csv")

################################################################################
# 1_6. Save metadata for visualization in Xenium Explorer
################################################################################

cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types){
final_result_select<- final_result[which(final_result$cancertype==cancer),]
final_result_select<-droplevels(final_result_select)
cell_types<-unique(final_result_select$celltype)
for (cell in cell_types){
final_result_select_cell<-final_result_select[which(final_result_select$celltype==cell),]
metadata_colonysize<-data.frame(cell_id=final_result_select_cell$Cell_ID,group=final_result_select_cell$size_group)
metadata_colonysize_colonyid<-data.frame(cell_id=final_result_select_cell$Cell_ID,group=final_result_select_cell$unique_colony_id)
write.csv(metadata_colonysize,paste0("~/space/analysis/colony/2_5_metadata_colonysize_",cancer,"_",cell,".csv"),row.names = FALSE)
write.csv(metadata_colonysize_colonyid,paste0("~/space/analysis/colony/2_5_metadata_colonysize_colonyid_",cancer,"_",cell,".csv"),row.names = FALSE)
}
}