library(Seurat)
library(dplyr)
library(purrr)
library(BiocParallel)
library(SeuratWrappers)
library(future)
library(data.table)
library(tidyverse)
options(future.globals.maxSize = 300 * 1024^3) 


#################1. Calculate the nearest distance between each malignant cell and T cell across cancer types
#################
#################
#################
#################
#################
#################


# Load metadata data
seurat_niches<-readRDS("~/space/seurat_niche/metadata_merge.rds") # Merged meta data of all cancers was generated in Step02, line 139

cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
# Define a container to store the results
result_df_merge <- tibble(
  Cell_ID = character(),
  cancertype = character(),
  sample = factor(),
  coord_x_malignant = numeric(),
  coord_y_malignant = numeric(),
  distance = numeric(),
  Cell_ID_tcell = character(),
  coord_x_tcell = numeric(),
  coord_y_tcell = numeric())

# Calculate the nearest distance between each malignant cell and T cell across cancer types
for (cancer in cancer_types){
  seurat_niches_select<-seurat_niches[which(seurat_niches$cancertype==cancer&seurat_niches$niche=="4"),]
  result_df <- seurat_niches_select %>%
    group_by(sample) %>%
    group_modify(~ {
      sample_data <- .x
      malignant_cells <- sample_data %>% filter(celltype == "Malignant")
      t_cells <- sample_data %>% filter(celltype == "T_cell")
      if(nrow(malignant_cells) == 0 | nrow(t_cells) == 0) {
        return(data.frame(
          Cell_ID = character(),
          cancertype = character(),
          coord_x_malignant = numeric(),
          coord_y_malignant = numeric(),
          distance = numeric(),
          Cell_ID_tcell = character(),
          coord_x_tcell = numeric(),
          coord_y_tcell = numeric()
        ))
      }
      # Loop for each cell
      result_list <- lapply(1:nrow(malignant_cells), function(i) {
        malignant_coords <- c(malignant_cells$coord_x[i], malignant_cells$coord_y[i])
        
        t_cell_distances <- sqrt(
          (t_cells$coord_x - malignant_coords[1])^2 + 
            (t_cells$coord_y - malignant_coords[2])^2
        )
        min_index <- which.min(t_cell_distances)
        # Return a list of results
        list(
          Cell_ID = as.character(malignant_cells$Cell_ID[i]),
          cancertype = as.character(malignant_cells$cancertype[i]),
          coord_x_malignant = as.numeric(malignant_cells$coord_x[i]),
          coord_y_malignant = as.numeric(malignant_cells$coord_y[i]),
          distance = as.numeric(t_cell_distances[min_index]),
          Cell_ID_tcell = as.character(t_cells$Cell_ID[min_index]),
          coord_x_tcell = as.numeric(t_cells$coord_x[min_index]),
          coord_y_tcell = as.numeric(t_cells$coord_y[min_index])
        )
      })
      # Combine the result list
      do.call(rbind.data.frame, result_list)
    }) %>%
    ungroup() %>%
    select(Cell_ID, cancertype, sample, coord_x_malignant, coord_y_malignant, 
           distance, Cell_ID_tcell, coord_x_tcell, coord_y_tcell)
  result_df_merge<-rbind(result_df_merge,result_df)
}
saveRDS(result_df_merge,"~/space/analysis/gradient/1_1_result_df_merge.rds")


################# 2. Identify gradient gene expression at the group level
#################
#################
#################
#################
#################
#################


################################################################################
# 2_1. Classify cells into different groups according to distance
################################################################################

max(result_df_merge$distance)
min(result_df_merge$distance)

result_df_grouped <- result_df_merge %>%
  group_by(cancertype) %>%
  arrange(distance, .by_group = TRUE) %>%  # Sort by distance
  mutate(group = ntile(distance, 10)) %>%   # Divide into 10 groups
  ungroup() %>%
  arrange(cancertype, group)  # Sort by cancer type and group
head(result_df_grouped, 10)
result_df_grouped$group<-paste0("group_",result_df_grouped$group)

saveRDS(result_df_grouped,"~/space/analysis/gradient/2_1_result_df_grouped.rds")

################################################################################
# 2_2. Calculate mean and median distance for each group in each cancer type
################################################################################

group_stats <- result_df_merge %>%
  group_by(cancertype) %>%
  arrange(distance, .by_group = TRUE) %>%
  mutate(group = ntile(distance, 10)) %>%
  group_by(cancertype, group) %>%
  summarise(
    n_cells = n(),
    mean_distance = mean(distance, na.rm = TRUE),
    median_distance = median(distance, na.rm = TRUE),
    sd_distance = sd(distance, na.rm = TRUE),
    min_distance = min(distance, na.rm = TRUE),
    max_distance = max(distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cancertype, group)
print(group_stats, n = Inf)
group_stats$group<-paste0("group_",group_stats$group)
saveRDS(group_stats,"~/space/analysis/gradient/2_2_group_stats.rds")

################################################################################
# 3_3. Load Seurat data for each cancer type
################################################################################

cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
# Load data of each cancer type
for(cancer in cancer_types){
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds"))
  assign(paste0("obj_", cancer), obj)
  message("Processed: ", cancer)
}

################################################################################
# 2_4. Identify distance gradient gene expression
################################################################################

library(WGCNA)
enableWGCNAThreads(nThreads = 8)
cor_df_grouped_merge<-data.frame()
for (cancer in cancer_types){
  obj<-get(paste0("obj_",cancer))
  result_df_select<-result_df_grouped[which(result_df_grouped$cancertype==cancer),]
  obj_select<-subset(obj,cells=result_df_select$Cell_ID)
  data_data <- GetAssayData(obj_select, assay = "Xenium", layer = "data")
  
  group_stats_select<-group_stats[which(group_stats$cancertype==cancer),]
  
  # Calculate the mean of each gene in each group
  common_cells <- intersect(result_df_select$Cell_ID, colnames(data_data))
  cell_order <- match(common_cells, result_df_select$Cell_ID)
  group_info <- result_df_select[cell_order, c("Cell_ID", "group")]
  data_filtered <- data_data[, common_cells]
  
  unique_groups <- unique(group_info$group)
  unique_groups <- sort(unique_groups)  # Ensure consistent order
  result_list <- lapply(unique_groups, function(g) {
    group_cells <- group_info$Cell_ID[group_info$group == g]
    valid_cells <- intersect(group_cells, colnames(data_filtered))
    if (length(valid_cells) > 0) {
      group_data <- data_filtered[, valid_cells, drop = FALSE]
      gene_means <- Matrix::rowMeans(group_data, na.rm = TRUE)
    } else {
      gene_means <- setNames(rep(NA, nrow(data_filtered)), rownames(data_filtered))
    }
    return(gene_means)
  })
  result_df_final <- as.data.frame(do.call(cbind, result_list))
  colnames(result_df_final) <- unique_groups
  
  # WGCNA analysis
  cor_results <- WGCNA::cor(
    x = group_stats_select$mean_distance,
    y = t(result_df_final[,as.character(group_stats_select$group)]),
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  # Convert to vector
  cor_vector <- as.numeric(cor_results)
  # Calculate sample size (effective sample size after considering missing values)
  n <- sum(!is.na(group_stats_select$mean_distance))  # Effective sample size of variable a
  # Calculate t-statistic and p-value
  t_stats <- cor_vector * sqrt((n - 2) / (1 - cor_vector^2))
  p_values <- 2 * pt(-abs(t_stats), df = n - 2)
  # Create result data frame
  cor_df_grouped <- data.frame(
    variable = rownames(result_df_final[,as.character(group_stats_select$group)]),
    correlation = cor_vector,
    t_statistic = t_stats,
    p_value = p_values,
    p_adjust = p.adjust(p_values, method = "fdr"),
    cancertype=rep(cancer,length(rownames(result_df_final[,as.character(group_stats_select$group)])))
  )
  cor_df_grouped_merge<-rbind(cor_df_grouped_merge,cor_df_grouped)
}
colnames(cor_df_grouped_merge)<-c("variable","correlation","t_statistic","p_value","p_adjust","cancertype")
write.csv(cor_df_grouped_merge,"~/space/analysis/gradient/2_4_cor_df_grouped_merge.csv")




################# 3. Identify gradient gene expression based on single-cell level
#################
#################
#################
#################
#################
#################


################################################################################
# 3_1. Identify distance gradient gene expression
################################################################################
library(WGCNA)
enableWGCNAThreads(nThreads = 8)

cor_df_merge<-data.frame()
for (cancer in cancer_types){
  # Select cancer type
  obj<-get(paste0("obj_",cancer))
  result_df_select<-result_df_merge[which(result_df_merge$cancertype==cancer),]
  obj_select<-subset(obj,cells=result_df_select$Cell_ID)
  data_data <- GetAssayData(obj_select, assay = "Xenium", layer = "data")
  # WGCNA analysis
  cor_results <- WGCNA::cor(
    x = result_df_select$distance,
    y = t(data_data[,as.character(result_df_select$Cell_ID)]),
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  # Convert to vector
  cor_vector <- as.numeric(cor_results)
  # Calculate sample size (effective sample size after considering missing values)
  n <- sum(!is.na(result_df_select$distance))  # Effective sample size of variable a
  # Calculate t-statistic and p-value
  t_stats <- cor_vector * sqrt((n - 2) / (1 - cor_vector^2))
  p_values <- 2 * pt(-abs(t_stats), df = n - 2)
  # Create result data frame
  cor_df <- data.frame(
    variable = rownames(data_data[,as.character(result_df_select$Cell_ID)]),
    correlation = cor_vector,
    t_statistic = t_stats,
    p_value = p_values,
    p_adjust = p.adjust(p_values, method = "fdr"),
    cancertype=rep(cancer,length(rownames(data_data[,as.character(result_df_select$Cell_ID)])))
  )
  cor_df_merge<-rbind(cor_df_merge,cor_df)
}
colnames(cor_df_merge)<-c("variable","correlation","t_statistic","p_value","p_adjust","cancertype")

# Svae the results
saveRDS(cor_df_merge,"~/space/analysis/gradient/3_1_cor_df_merge.rds")
write.csv(cor_df_merge,"~/space/analysis/gradient/3_1_cor_df_merge.csv")