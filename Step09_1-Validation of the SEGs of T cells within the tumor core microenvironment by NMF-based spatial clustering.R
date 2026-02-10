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

################# 1. Identify neighboring cell and calculate neighboring cell composition in each cancer type
#################
#################
#################
#################
#################
#################

################################################################################
# 1_1. Configuration section
################################################################################

# Set parameters
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")

# Input parameters
input_dir <- "~/space/metadata"
output_dir <- "~/space/analysis/neighbor/allcancer/tcell"

# Target cell type for neighborhood analysis
target_celltype <- "T_cell"

# Distance threshold for neighbor identification (in microns)
neighbor_radius <- 80

# Parallel processing settingsDefine of the main analysis function
n_cores <- 16

################################################################################
# 1_2. Load required libraries
################################################################################
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(doParallel)
})

################################################################################
# 1_3. Define the main analysis function
################################################################################
analyze_cancer_neighborhood <- function(cancer_type, input_dir, output_dir, 
                                        target_celltype = "T_cell", 
                                        neighbor_radius = 80, 
                                        n_cores = 16) {
  
  cat("\n" + strrep("=", 80) + "\n")
  cat("Starting analysis for:", toupper(cancer_type), "\n")
  start_time <- Sys.time()
  cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("-", 80) + "\n")
  
  # (1). Load and prepare the metadata
  cat("1. Loading metadata...\n")
  input_file <- file.path(input_dir, paste0("obj_", cancer_type, "_metadata.rds")) # Meta data of each cancer were generated in Step01, line 376
  
  if (!file.exists(input_file)) {
    cat("Warning: Input file not found:", input_file, "\n")
    return(NULL)
  }
  
  # Read data and convert to data.table
  metadata <- as.data.table(readRDS(input_file), keep.rownames = "cellid")
  
  cat("  Total cells:", nrow(metadata), "\n")
  cat("  Cell types:", paste(sort(unique(metadata$celltype)), collapse = ", "), "\n")
  
  # (2). Filter target cells
  cat("2. Filtering target cells (", target_celltype, ")...\n", sep = "")
  metadata_target <- metadata[celltype == target_celltype, ]
  
  if (nrow(metadata_target) == 0) {
    cat("Warning: No target cells found for", cancer_type, "\n")
    return(NULL)
  }
  
  cat("  Target cells:", nrow(metadata_target), "\n")
  cat("  Number of samples:", length(unique(metadata_target$sample)), "\n")
  
  # Get all cell types in dataset
  all_celltypes <- sort(unique(metadata$celltype))
  
  # (3). Setup parallel processing
  cat("3. Setting up parallel processing (", n_cores, " cores)...\n", sep = "")
  registerDoParallel(cores = n_cores)
  
  # (4). Process by sample
  cat("4. Identifying neighbors for each target cell...\n")
  samples <- unique(metadata_target$sample)
  cat("  Processing", length(samples), "samples\n")
  
  results <- foreach(i = samples, .combine = rbind, .packages = "data.table") %dopar% {
    # Get data for current sample
    sample_target <- metadata_target[sample == i, ]
    sample_all <- metadata[sample == i, ]
    coords_all <- as.matrix(sample_all[, .(coord_x, coord_y)])
    
    # Process each target cell
    sample_neib <- rbindlist(lapply(1:nrow(sample_target), function(j) {
      # Get target cell coordinates
      target_coords <- as.numeric(sample_target[j, .(coord_x, coord_y)])
      
      # Calculate Euclidean distances (vectorized)
      dists <- sqrt((coords_all[,1] - target_coords[1])^2 + 
                      (coords_all[,2] - target_coords[2])^2)
      
      # Identify neighbors within radius (excluding self)
      neighbor_idx <- which(dists <= neighbor_radius & dists > 0)
      
      if (length(neighbor_idx) == 0) return(NULL)
      
      # Extract neighbor information
      neighbors <- sample_all[neighbor_idx, ]
      neighbors[, `:=`(
        dist = dists[neighbor_idx], 
        targetcellid = sample_target[j, cellid]
      )]
      
      return(neighbors)
    }), fill = TRUE)
    
    return(sample_neib)
  }
  
  cat("  Total neighbor pairs identified:", nrow(results), "\n")
  
  # (5). Calculate the cell type composition
  cat("5. Calculating cell type composition...\n")
  
  if (nrow(results) > 0) {
    # Count neighbors by cell type for each target cell
    cellnumber <- dcast(
      results[, .N, by = .(celltype, targetcellid)],
      celltype ~ targetcellid, 
      value.var = "N", 
      fill = 0
    )
    
    # Ensure all cell types are represented
    missing_types <- setdiff(all_celltypes, cellnumber$celltype)
    if (length(missing_types) > 0) {
      # Create zero rows for missing cell types
      zero_matrix <- matrix(0, nrow = length(missing_types), 
                            ncol = ncol(cellnumber) - 1)
      colnames(zero_matrix) <- names(cellnumber)[-1]
      
      cellnumber <- rbind(
        cellnumber,
        data.table(celltype = missing_types, zero_matrix)
      )
    }
    
    # Calculate proportions
    cellcomposition <- as.data.frame(cellnumber[, -1, with = FALSE])
    
    if (ncol(cellcomposition) > 0) {
      cellcomposition <- sweep(cellcomposition, 2, colSums(cellcomposition), "/")
    }
    
    rownames(cellcomposition) <- cellnumber$celltype
    
  } else {
    cat("  No neighbors found. Creating empty results.\n")
    # Create empty results structure
    cellnumber <- data.table(celltype = all_celltypes, 
                             matrix(0, nrow = length(all_celltypes), ncol = 0))
    cellcomposition <- data.frame(row.names = all_celltypes)
  }
  
  # (6). Save results
  cat("6. Saving results...\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save RDS files
  saveRDS(cellnumber, file.path(output_dir, paste0("cellnumber_", cancer_type, ".rds")))
  saveRDS(results, file.path(output_dir, paste0("data_neib_", cancer_type, ".rds")))
  saveRDS(cellcomposition, file.path(output_dir, paste0("cellcomposition_", cancer_type, ".rds")))
  
  # Save CSV file
  if (nrow(cellcomposition) > 0 && ncol(cellcomposition) > 0) {
    write.csv(cellcomposition, file.path(output_dir, paste0("cellcomposition_", cancer_type, ".csv")))
  }
  
  # (7). Generate summary statistics
  end_time <- Sys.time()
  cat("\n" + strrep("-", 80) + "\n")
  cat("Summary for", toupper(cancer_type), ":\n")
  cat("  Total target cells:", nrow(metadata_target), "\n")
  cat("  Total neighbor pairs:", nrow(results), "\n")
  if (nrow(results) > 0) {
    cat("  Average neighbors per target cell:", 
        round(nrow(results)/nrow(metadata_target), 2), "\n")
  }
  cat("  Unique cell types in neighborhood:", length(unique(results$celltype)), "\n")
  cat("  Analysis completed:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("  Total time:", format(end_time - start_time, digits = 2), "\n")
  cat(strrep("=", 80) + "\n")
  
  # Return summary results
  return(list(
    cancer_type = cancer_type,
    n_target_cells = nrow(metadata_target),
    n_neighbor_pairs = nrow(results),
    n_cell_types = length(all_celltypes),
    processing_time = end_time - start_time
  ))
}

################################################################################
# 1_4. Main execution loop
################################################################################
main <- function() {
  cat(strrep("#", 80) + "\n")
  cat("SPATIAL NEIGHBORHOOD ANALYSIS PIPELINE\n")
  cat("Target cell type:", target_celltype, "\n")
  cat("Neighbor radius:", neighbor_radius, "microns\n")
  cat("Number of cores:", n_cores, "\n")
  cat("Number of cancer types:", length(cancer_types), "\n")
  cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("#", 80) + "\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Initialize a results lists
  all_results <- list()
  successful_runs <- 0
  failed_runs <- 0
  
  # Process for each cancer type
  for (cancer in cancer_types) {
    tryCatch({
      result <- analyze_cancer_neighborhood(
        cancer_type = cancer,
        input_dir = input_dir,
        output_dir = output_dir,
        target_celltype = target_celltype,
        neighbor_radius = neighbor_radius,
        n_cores = n_cores
      )
      
      if (!is.null(result)) {
        all_results[[cancer]] <- result
        successful_runs <- successful_runs + 1
      } else {
        failed_runs <- failed_runs + 1
        cat("Analysis failed for:", cancer, "\n")
      }
      
    }, error = function(e) {
      cat("\nERROR processing", cancer, ":\n")
      cat("  ", conditionMessage(e), "\n")
      failed_runs <- failed_runs + 1
    })
  }

################################################################################
# 1_5. Generate summary report
################################################################################
  cat("\n" + strrep("#", 80) + "\n")
  cat("ANALYSIS COMPLETE\n")
  cat(strrep("-", 80) + "\n")
  cat("Summary:\n")
  cat("  Successful runs:", successful_runs, "/", length(cancer_types), "\n")
  cat("  Failed runs:", failed_runs, "/", length(cancer_types), "\n")
  cat("  End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  if (length(all_results) > 0) {
    cat("\nDetailed summary by cancer type:\n")
    cat(strrep("-", 40) + "\n")
    for (result in all_results) {
      cat(sprintf("%-8s: %4d target cells, %6d neighbor pairs, %3d cell types, %6.1f seconds\n",
                  toupper(result$cancer_type),
                  result$n_target_cells,
                  result$n_neighbor_pairs,
                  result$n_cell_types,
                  as.numeric(result$processing_time, units = "secs")))
    }
  }
  
  # Save summary statistics
  summary_file <- file.path(output_dir, "analysis_summary.txt")
  sink(summary_file)
  cat("Spatial Neighborhood Analysis Summary\n")
  cat("====================================\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Target cell type:", target_celltype, "\n")
  cat("Neighbor radius:", neighbor_radius, "microns\n")
  cat("\nResults by cancer type:\n")
  cat(strrep("-", 60), "\n")
  for (result in all_results) {
    cat(sprintf("%-8s: Targets=%4d, Neighbors=%6d, CellTypes=%2d, Time=%5.1fs\n",
                toupper(result$cancer_type),
                result$n_target_cells,
                result$n_neighbor_pairs,
                result$n_cell_types,
                as.numeric(result$processing_time, units = "secs")))
  }
  sink()
  
  cat("\nSummary saved to:", summary_file, "\n")
  cat(strrep("#", 80) + "\n")
  
  return(all_results)
}

################################################################################
# 1_6. Execute analysis
################################################################################

library(data.table)
library(parallel)
library(doParallel)
registerDoParallel(cores = 16)
base_metadata_dir <- "~/space/metadata/"
base_result_dir <- "~/space/analysis/neighbor/allcancer/tcell/"  # Result saving directory


for (cancer_type in cancer_types) {

  # (1). Print current processing progress and start time

  cat("=====================================\n")
  cat("Processing cancer type:", cancer_type, "\n")
  start_time <- Sys.time()
  cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  

  # (2). Construct dynamic file path and check if the file exist

  metadata_file <- paste0(base_metadata_dir, "obj_", cancer_type, "_metadata.rds") # Meta data of each cancer were generated in Step01, line 376
  if (!file.exists(metadata_file)) {
    message("Error: Metadata file for ", cancer_type, " does not exist. Skipping...\n")
    next
  }
  

  # (3). Read and preprocess metadata

  # Convert to data.table format
  metadata <- as.data.table(readRDS(metadata_file), keep.rownames = "cellid")
  
  # Filter target cell type (T_cell)
  metadata_target <- metadata[celltype == "T_cell", ]
  all_celltypes <- sort(unique(metadata$celltype))
  

  # (4). Process sample group data in parallel

  samples <- unique(metadata_target$sample)
  results <- foreach(i = samples, .combine = rbind, .packages = "data.table") %dopar% {
    # Get data for the current sample
    sample_target <- metadata_target[sample == i, ]
    sample_all <- metadata[sample == i, ]
    coords_all <- as.matrix(sample_all[, .(coord_x, coord_y)])
    
    # Process each target cell
    sample_neib <- rbindlist(lapply(1:nrow(sample_target), function(j) {
      target_coords <- as.numeric(sample_target[j, .(coord_x, coord_y)])
      
      # Vectorized distance calculation
      dists <- sqrt((coords_all[,1] - target_coords[1])^2 + 
                      (coords_all[,2] - target_coords[2])^2)
      
      # Define neighboring cells (Exclude self, distance ≤ 80)
      neighbor_idx <- which(dists <= 80 & dists > 0)
      if (length(neighbor_idx) == 0) return(NULL)
      
      # Return neighboring cell information
      neighbors <- sample_all[neighbor_idx, ]
      neighbors[, `:=`(dist = dists[neighbor_idx], 
                       targetcellid = sample_target[j, cellid])]
      return(neighbors)
    }))
    
    return(sample_neib)
  }
  

  # (5). Subsequent data processing

  # Count the number of neighboring cell types
  cellnumber <- dcast(results[, .N, by = .(celltype, targetcellid)],
                      celltype ~ targetcellid, value.var = "N", fill = 0)
  
  # Ensure all cell types are included
  missing_types <- setdiff(all_celltypes, cellnumber$celltype)
  if (length(missing_types) > 0) {
    cellnumber <- rbind(cellnumber,
                        data.table(celltype = missing_types,
                                   matrix(0, nrow = length(missing_types),
                                          ncol = ncol(cellnumber) - 1)))
  }
  
  # Calculate cell composition proportions
  cellcomposition <- as.data.frame(cellnumber[, -1, with = FALSE])
  cellcomposition <- sweep(cellcomposition, 2, colSums(cellcomposition), "/")
  rownames(cellcomposition) <- cellnumber$celltype
  

  # (6). Save results

  # Construct paths for all result files
  cellnumber_file <- paste0(base_result_dir, "cellnumber_", cancer_type, ".rds")
  neib_data_file <- paste0(base_result_dir, "data_neib_", cancer_type, ".rds")
  cellcomp_rds_file <- paste0(base_result_dir, "cellcomposition_", cancer_type, ".rds")
  cellcomp_csv_file <- paste0(base_result_dir, "cellcomposition_", cancer_type, ".csv")
  
  # Save results in rds and csv format
  saveRDS(cellnumber, cellnumber_file)
  saveRDS(results, neib_data_file)
  saveRDS(cellcomposition, cellcomp_rds_file)
  write.csv(cellcomposition, cellcomp_csv_file, row.names = TRUE)
  

  end_time <- Sys.time()
  cat("Processing completed for", cancer_type, "\n")
  cat("Script end time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total running time for", cancer_type, ":", format(end_time - start_time, digits = 2), "\n\n")
}

stopImplicitCluster()  # Shut down the parallel computing cluster and release resources
cat("All cancer types processed! All results saved to:", base_result_dir, "\n")


################# 2. Merge cell composition data
#################
#################
#################
#################
#################
#################
library(data.table)
library(dplyr)
library(tibble)
library(purrr)

# (1). Set cancer type names and input directory
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
base_path <- "~/space/analysis/neighbor/allcancer/tcell/"

# (2). Extract data
cellcomposition_list <- lapply(cancer_types, function(type) {
  file_path <- paste0(base_path, "cellcomposition_", type, ".csv")
  dt <- fread(file_path, encoding = "UTF-8")  # fread does not retain row names by default, manual processing required
  df <- as.data.frame(dt)
  rownames(df) <- df[, 1]  # Assume the first column is row names
  df <- df[, -1, drop = FALSE]# Remove the first column (row names already stored)
  return(df)
})
names(cellcomposition_list) <- cancer_types
length(cellcomposition_list)
cellcomposition_list[["paad"]][,1:6]

# (3). Check data
lapply(cellcomposition_list, dim)  # Check dimensions
lapply(cellcomposition_list, function(df) df[1:5, 1:5])  # View first 5 rows and 5 columns

# (4). Merge cell composition data
long_list <- lapply(cancer_types, function(type) {
  rownames_to_column(as.data.frame(cellcomposition_list[[type]]), var = "CellType")
})
names(long_list) <- paste0(cancer_types, "_long")
list_of_dfs <- long_list

# (5). Use reduce() to progressively full_join by "CellType"
cellcomposition_merge <- reduce(list_of_dfs, full_join, by = "CellType")
dim(cellcomposition_merge)

# (6). Set Cell-type back as row names
cellcomposition_merge<-as.data.frame(cellcomposition_merge)
rownames(cellcomposition_merge) <- cellcomposition_merge$CellType
cellcomposition_merge <- cellcomposition_merge %>% select(-CellType)

# (7). Replace NA with 0
setnafill(cellcomposition_merge, fill = 0)
dim(cellcomposition_merge)
ncol(cellcomposition_merge)==(sum(sapply(list_of_dfs, ncol))-length(cancer_types))

# (8). Save data
cellcomposition_merge[1:5,1:5]
cellcomposition_merge=cellcomposition_merge[c("B_cell","Endothelial","Fibroblast","Malignant","Myeloid","Smooth_muscle","Neural","Pericyte","Plasma","T_cell","Oligodendrocyte","Keratinocyte"),]
write.csv(cellcomposition_merge,"~/space/analysis/neighbor/allcancer/tcell/cellcomposition_merge.csv")
saveRDS(cellcomposition_merge,"~/space/analysis/neighbor/allcancer/tcell/cellcomposition_merge.rds")
