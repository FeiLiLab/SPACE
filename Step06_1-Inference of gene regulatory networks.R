library(SeuratDisk)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(data.table)    
library(SummarizedExperiment)
options(future.globals.maxSize = 300 * 1024^3)



################# 1. Prepare data for pySCENIC transcription factor (TF) analysis
#################
#################
#################
#################
#################
#################

################################################################################
# 1_1. Set parameters
################################################################################

# Set work directory
setwd("~/space")
dir.create("pyscenic")
setwd("~/space/pyscenic")
output_dir <- "./loom_files"
dir.create(output_dir, showWarnings = FALSE)

# Define cell types to exclude
cells_to_remove <- c("Malignant", "Oligodendrocyte", "Keratinocyte")

################################################################################
# 1_2. Convert the Seurat object to loom file for pySCENIC input
################################################################################
cancer_list <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")

for (cancer in cancer_list) {

  # Read the RDS file
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS files of each cancer were generated in Step02, line 135
  
  # Filter out unwanted cell types
  obj_filtered <- subset(obj, subset = !(celltype %in% cells_to_remove))
  
  # Randomly subset 1/3 of the remaining cells
  set.seed(123)
  n_cells <- ncol(obj_filtered)
  n_subset <- round(n_cells / 3)
  
  # Ensure there are cells left after filtering
  if (n_subset > 0) {
    selected_cells <- sample(Cells(obj_filtered), size = n_subset)
    obj_subset <- subset(obj_filtered, cells = selected_cells)
    
    # Define output path
    output_path <- paste0("~/space/pyscenic/loom_files/obj_",cancer, "_subset_filter.loom")
    
    # Skip failed files via tryCatch
    tryCatch({
      SaveLoom(
        object = obj_subset,
        filename = output_path,
        overwrite = FALSE
      )
      cat("Successfully saved:", output_path, "\n\n")
    }, error = function(e) {
      cat("Error saving", dataset_id, ":", conditionMessage(e), "\n")
    })
    
  } else {
    cat("Warning: No cells left for", dataset_id, "after filtering.\n\n")
  }

}


################# 2. Execute pySCENIC TF analysis via Linux Slurm file
#################
#################
#################
#################
#################
#################

# Notice: This step was performed with Liunx Slurm file of pySCENIC analysis: Step06_2-Inference of gene regulatory networks.slurm

################# 3. Pre-process TF analysis results for downstream visualization
#################
#################
#################
#################
#################
#################

################################################################################
# 3_1. Load results of pySCENIC TF analysis
################################################################################

# (1). Set parameters
setwd("~/space/pyscenic/analysis/")
cancer_list <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")

# (2). Initialize containers for pySCENIC output objects
loom_list       <- list()      
regulons_list      <- list()      
regulonAUC_list    <- list()      

# (3). Load results of each cancer
for (i in cancer_list) {
  # Set the file path for the pySCENIC output files
  loom_path <- paste0("~/space/pyscenic/analysis/obj_", i, "/obj_",i,"_aucell.loom") # The pySCENIC output files were generated with Step06_2-Inference of gene regulatory networks.slurm
  
  # Check if the loom file exists before attempting to open
  if (!file.exists(loom_path)) {
    warning(paste("Loom file not found for", i, ":", loom_path))
    next
  }
  cat("Processing:", i, "\n")
  
  # Open the loom connection
  loom <- open_loom(loom_path)
  loom_list[[i]] <- loom  
  
  # Extract the regulon incidence matrix (binary matrix: cell x regulon)
  # pySCENIC loom files usually use the attribute name "Regulons"
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  cat("Regulons matrix dim:", dim(regulons_incidMat), "\n")
  print(regulons_incidMat[1:4, 1:4])
  
  # Convert the incidence matrix to a gene list (target genes for each regulon)
  regulons_list[[i]] <- regulonsToGeneLists(regulons_incidMat)
  
  # Extract the AUC values
  regulonAUC_list[[i]] <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  
  # Close the loom connection to release the file lock
  close_loom(loom)
}


################################################################################
# 3_2. Construct Seurat objects from regulon activity matrices
################################################################################

# (1). Initialize empty lists for analysis results
seurat_list <- list()

# (2). Use a loop to create Seurat objects from the input data files
for (tumor in names(regulonAUC_list)) {
  
  cat("Processing tumor:", tumor, "\n")
  auc_mat <- regulonAUC_list[[tumor]]
  
  # Construct Seurat assay from regulon activity matrices
  auc_mat <- assay(auc_mat)
  
  # Identify cells present in both the AUC matrix and the metadata
  common_cells <- intersect(colnames(auc_mat), rownames(allcancer_combined_meta))
  
  if (length(common_cells) == 0) {
    stop(paste("No overlapping cells between AUC matrix and metadata for tumor", tumor))
  }
  
  cat("  Found", ncol(auc_mat), "regulons and", length(common_cells), "cells in common\n")
  
  # Subset the AUC matrix to keep only cells with available metadata
  auc_mat <- auc_mat[, common_cells, drop = FALSE]
  
  # Create a Seurat assay object
  auc_assay <- CreateAssayObject(data = auc_mat) 
  
  # Initialize the Seurat object using the AUC assay
  seu <- CreateSeuratObject(counts = auc_assay, assay = "AUC", meta.data = NULL)
  
  # Subset the global metadata for the current cancer
  seu_meta <- allcancer_combined_meta[common_cells, , drop = FALSE]
  
  # Add source information to the metadata for tracking
  seu_meta$orig.tumor <- tumor
  seu_meta$orig.ident <- tumor
  
  # Incorporate the metadata into the Seurat object
  seu <- AddMetaData(seu, metadata = seu_meta)
  
  # Append the processed object to the list
  seurat_list[[tumor]] <- seu
}

# (3) Save the result
saveRDS(seurat_list, file = "~/space/pyscenic/analysis/TF_RegulonAUC.rds")


################################################################################
# 3_3. Load metadata
################################################################################

# (1). Load metadata
allcancer_combined_meta=readRDS("~/space/seurat_niche/metadata_merge.rds") # Merged meta data of all cancers was generated in Step02, line 139
table(allcancer_combined_meta$majorcelltype)

# (2). Subset the metadata
all_pyscenic_cells=sapply(seurat_list, function(x){c(colnames(x))})
all_pyscenic_cells=unlist(all_pyscenic_cells)
allcancer_combined_meta_sub=allcancer_combined_meta[row.names(allcancer_combined_meta) %in% all_pyscenic_cells,]
allcancer_combined_meta_sub$cancer=tolower(allcancer_combined_meta_sub$cancer)


################# 4. Identify spatially enriched TFs (SETFs) and cross-cancer shared SETFs
#################
#################
#################
#################
#################
#################

################################################################################
# 4_1. Identify spatially enriched TFs
################################################################################

# (1). Set parameters
plan("multisession", workers = 16)
setf_all_target <- data.frame()
# Exclude cancer specific cell lineages
majorcelltype_list <- setdiff(names(table(allcancer_combined_meta_sub$majorcelltype)), c("Keratinocyte", "Oligodendrocyte"))

# (2). Identify spatially enriched TFs
for (majorcelltype_select in majorcelltype_list) {
  # Subset metadata for the currently selected cell type
  allcancer_combined_meta_sub_target <- allcancer_combined_meta_sub[which(allcancer_combined_meta_sub$majorcelltype == majorcelltype_select), ]
  
  # Filter cancer types where each niche has more than 2 cells to ensure statistical validity
  cancer_types <- allcancer_combined_meta_sub_target %>%
    group_by(cancer, niches) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(cancer) %>%
    summarise(all_above_3 = all(count > 2)) %>%
    filter(all_above_3) %>%
    pull(cancer)
  
  # Iterate through each valid cancer type to identify spatially enriched TFs
  for (cancer in cancer_types) {
    obj <- seurat_list[[cancer]]
    
    niche_list <- unique(obj$niches)
    Idents(obj) <- obj$niches
    
    for (niche_select in niche_list) {
      # Identify candidate SETFs for the current niche versus all other niches in the same cell type
      setf_all_target_1 <- FindMarkers(obj, ident.1 = niche_select, ident.2 = setdiff(niche_list, niche_select))
      
      # Append metadata information to the results table
      setf_all_target_1$cancer <- rep(cancer, nrow(setf_all_target_1))
      setf_all_target_1$niches <- rep(niche_select, nrow(setf_all_target_1))
      setf_all_target_1$gene <- row.names(setf_all_target_1)
      setf_all_target_1$majorcelltype <- rep(majorcelltype_select, nrow(setf_all_target_1))
      
      # Combine specific niche results into the data frame
      setf_all_target <- rbind(setf_all_target, setf_all_target_1)
    }
  }
}

# (3). Filter for up-regulated TFs
setf_all_target <- setf_all_target[setf_all_target$avg_log2FC > 0, ]
setf_all_target <- setf_all_target[setf_all_target$p_val_adj < 0.05, ]


################################################################################
# 4_2. Identify cancer-shared SETFs
################################################################################

# (1). Set parameters
majorcelltype_list <- names(table(setf_all_target$majorcelltype))
df_crosscancer_merge <- data.frame()

# (2). Identify cancer-shared SETFs
for (majorcelltype_select in majorcelltype_list) {
  # Filter data for the specific cell type
  setf_final_select <- setf_all_target[which(setf_all_target$majorcelltype == majorcelltype_select), ]
  
  # Identify all unique niches within this cell type
  niche_list <- names(table(setf_final_select$niche))
  
  for (niche_select in niche_list) {
    # Filter for the combination of niche and cell type
    setf_final_select_sub <- setf_all_target[which(setf_all_target$niche == niche_select & setf_all_target$majorcelltype == majorcelltype_select), ]
    
    # Calculate frequency of each TF across cancer types
    df_crosscancer <- as.data.frame(table(setf_final_select_sub$gene)[order(table(setf_final_select_sub$gene), decreasing = TRUE)])
    
    # Append niche and cell type metadata to the frequency table
    df_crosscancer$niche <- rep(niche_select, nrow(df_crosscancer))
    df_crosscancer$majorcelltype <- rep(majorcelltype_select, nrow(df_crosscancer))
    
    # Merge the current subset results into the master frequency data frame
    df_crosscancer_merge <- rbind(df_crosscancer_merge, df_crosscancer)
  }
}

################################################################################
# 4_3. Retain SETFs that were consistently detected (Percentage > 5%) in the corresponding cell type in the public pan-cancer scRNA-seq dataset
################################################################################

# (1). Remove the "(+)" suffix from regulon names to match standard gene symbols
setf_all_target$gene <- gsub("\\(\\+\\)", "", setf_all_target$gene)

# (2). Load the pan-cancer cell lineage gene expression percentage data
# Gene expression percentage data were derived from the public pan-cancer scRNA-seq dataset reported in "Cross-tissue multicellular coordination and its rewiring in cancer"
pancancer_gene <- read.csv("~/space/pancancer_celllineage.genes.percent.csv")

# (3). Identify overlapping cell types between the expression reference and the TF results
celltype_overlap <- intersect(names(table(setf_all_target$majorcelltype)), names(table(pancancer_gene$majorcelltype)))

# (4). Exclude specific myeloid lineages from the analysis
celltype_analysis <- setdiff(celltype_overlap, c("Monocyte", "Neutrophil"))

# (5). Initialize an empty data frame for the filtered results
setf_final <- data.frame()

# (6). Filter the candidate SETFs
for (celltype in celltype_analysis) {
  # Filter for genes expressed in more than 5% of the cells in the reference
  pancancer_gene_select <- pancancer_gene[which(pancancer_gene$majorcelltype == celltype & pancancer_gene$Proportion > 0.05), ]
  gene_expressed <- pancancer_gene_select$Gene 
  
  # Subset the results for the current cell type, keeping only validly expressed genes
  setf_final_1 <- setf_all_target[which(setf_all_target$majorcelltype == celltype & setf_all_target$gene %in% gene_expressed), ]
  
  # Combine the filtered subset into the master data frame
  setf_final <- rbind(setf_final, setf_final_1)
}

# (7). Export the final filtered results to a CSV file
write.csv(setf_final, "~/space/pyscenic/analysis/4_3_setf_all_target_ref0.05.csv")


################################################################################
# 4_4. Retain cancer-shared SETFs that were consistently detected (Percentage > 5%) in the corresponding cell type in the public pan-cancer scRNA-seq dataset
################################################################################

# (1). Clean regulon names by removing the (+) suffix to get standard gene symbols
df_crosscancer_merge$gene <- gsub("\\(\\+\\)", "", df_crosscancer_merge$Var1)

# (2). Identify overlapping cell types between the expression reference and the TF results
celltype_overlap <- intersect(names(table(df_crosscancer_merge$majorcelltype)), names(table(pancancer_gene$majorcelltype)))

# (3). Exclude specific myeloid lineages from the analysis
celltype_analysis <- setdiff(celltype_overlap, c("Monocyte", "Neutrophil"))

# (4). Initialize an empty data frame for the final filtered results
setf_final_crosscancer <- data.frame()

# (5). Filter the candidate cancer-shared SETFs
for (celltype in celltype_analysis) {
  # Filter for genes expressed in more than 5% of the cells in the reference
  pancancer_gene_select <- pancancer_gene[which(pancancer_gene$majorcelltype == celltype & pancancer_gene$Proportion > 0.05), ]
  gene_expressed <- pancancer_gene_select$Gene 
  
  # Filter cancer-shared results to keep only those genes that meet the expression threshold
  setf_final_1 <- df_crosscancer_merge[which(df_crosscancer_merge$majorcelltype == celltype & df_crosscancer_merge$gene %in% gene_expressed), ]
  
  # Append the filtered results into the data frame
  setf_final_crosscancer <- rbind(setf_final_crosscancer, setf_final_1)
}

# (6). Save the cancer-shared TF frequency table
write.csv(setf_final_crosscancer, "~/space/pyscenic/analysis/4_4_setf_crosscancer_merge_ref0.05.csv")


