options(future.globals.maxSize = 300 * 1024^3)
library(Seurat)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(Cairo)
library(arrow)
library(ggrepel)
library(future)
library(data.table)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(tidyr)

colorpanel <- c("Fibroblast" = '#A6CEE3',"Neural" = '#1F78B4',"T_cell" = '#B53E2B',"Myeloid" = '#33A02C',
                "Myofibroblast" = '#FB9A99',"Plasma" = '#B2DF8A',"Endothelial" = '#FF7F00',"Pericyte" = '#6A3D9A',
                "Malignant/Epithelial" = '#CAB2D6',"B_cell" = '#FDBF6F', 
                "Oligodendrocyte" = "#FFFF99","Keratinocyte"="#FEE1D2")

colorniche <- c("1" = '#466791',"2" = '#E39A35',"3" = '#712820',"4" = '#C1E6F3',
                "5" = '#df462a', "6" = "#CAB2D6", "7" = "#60bf37", "8" = '#6778AE')

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
mypal_1 <- union(pal_npg("nrc", alpha = 0.7)(10),my36colors)
mypal_2<-c('#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#B53E2B','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928'
           ,"#BBBDDC","#C1E6F3","#FEE1D2","#7DADC6","#F7A073","#E6754B","#AD7F7B")
mypal_3<-union(mypal_2,mypal_1)
mypal_4<-unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))[-c(6,17,18,19)][-c(7)]
mypal_5<-c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
mypal_6<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

color_cancertype <- c("BLCA" = "#8DD3C7", "BRCA" = "#FFFFB3", "CRC" = "#BEBADA", "ESCC" = "#FB8072", 
                      "GBM" = "#80B1D3", "LIHC" = "#FDB462", "NSCLC" = "#B3DE69", "PAAD" = "#FCCDE5", 
                      "PRAD" = "#D9D9D9", "RCC" = "#BC80BD", "SKCM" = "#CCEBC5", "STAD" = "#FFED6F")



################# 1. Correlation analysis of distance-based gradient gene expression
#################
#################
#################
#################
#################
#################

################################################################################
# 1_1. Prepare data for gradient gene analysis
################################################################################

# (1). Set samples to perform gradient analysis 
gradient_caners=c("BLCA_04","BLCA_05","BLCA_07","BLCA_12","BLCA_17","BLCA_29","BLCA_30","BLCA_31",
                  "BLCA_34","BLCA_37","BLCA_39","BRCA_03","BRCA_10","BRCA_11","BRCA_13","BRCA_16","BRCA_17","BRCA_20","BRCA_26",
                  "BRCA_29","BRCA_31","BRCA_37","CRC_02","CRC_06","CRC_08","CRC_09","CRC_12","CRC_17","CRC_22","CRC_23",
                  "CRC_24","CRC_25","CRC_28","CRC_34","CRC_37","ESCC_01","ESCC_02","ESCC_07","ESCC_09","ESCC_10","ESCC_11","ESCC_12","ESCC_22",
                  "ESCC_27","ESCC_31","ESCC_32","ESCC_36","ESCC_37","GBM_19","GBM_30","GBM_31","LIHC_07","LIHC_08","LIHC_11","LIHC_13","LIHC_14","LIHC_15","LIHC_16","LIHC_17",
                  "LIHC_20","LIHC_30","LIHC_32","NSCLC_03","NSCLC_04","NSCLC_05","NSCLC_06","NSCLC_07","NSCLC_13",
                  "NSCLC_16","NSCLC_17","PAAD_01","PAAD_03","PAAD_12","PAAD_24","PAAD_26","PAAD_34","PRAD_02","PRAD_18","PRAD_24","PRAD_25","PRAD_26","PRAD_29",
                  "PRAD_31","PRAD_38","PRAD_42","RCC_02","RCC_05","RCC_11","RCC_12","RCC_13","RCC_17","RCC_19","RCC_24","RCC_28",
                  "RCC_29","RCC_38","RCC_40","SKCM_03","SKCM_08","SKCM_09","SKCM_12","SKCM_14","SKCM_17","SKCM_19","SKCM_20","SKCM_21","SKCM_24","SKCM_26","SKCM_28",
                  "SKCM_29","SKCM_30","SKCM_31","SKCM_36","SKCM_39","SKCM_41","STAD_09","STAD_15","STAD_17","STAD_20","STAD_21","STAD_23","STAD_26","STAD_31","STAD_35",
                  "STAD_37","STAD_39")

# (2). Load data from RDS
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types) {
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS files of each cancer were generated in Step02, line 135
  colnames(obj@assays$Xenium) = row.names(obj@meta.data)
  Idents(obj) <- obj$celltype
  obj@images = list()
  gradient_idx <- intersect(gradient_caners, unique(obj$sample))
  obj <- subset(obj, sample %in% gradient_idx)
  assign(paste0("obj_",cancer,"_gradient_caners"),obj)
}


# (3). Load boundary data
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC", "GBM", "LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM", "STAD")
dist_list <- list()
for (ct in cancer_types) {
  base_path <- paste0("~/space/boundary/boundary_analysis_results", ct, "/boundary_with_zones.csv")
  
  message(">>> Processing: ", ct)
  df <- read.csv(base_path)
  
  # The distance on the tumor side is set to negative to distinguish between the inside and outside of the tumor
  df$dist_to_boundary[df$zoning_longest == "Epithelial_Zone"] <- 
    -(df$dist_to_boundary[df$zoning_longest == "Epithelial_Zone"])
  
  # Match cell id
  lower_ct <- tolower(ct)
  df$cellid <- gsub(paste0("_", lower_ct), paste0("_", ct), df$cellid)
  df$sample <- gsub(paste0(lower_ct, "_"), paste0(ct, "_"), df$sample)
  
  dist_list[[ct]] <- df
}


# (4). Merge boundary and distance data to metadata and combine metadata of all cancer types
all_meta_list <- list()

for (ct in cancer_types) {
  message(">>> Processing Metadata for: ", ct)
  
  # Current meta data
  obj_name <- paste0("obj_", ct, "_gradient_caners")
  curr_obj <- get(obj_name)

  # Current boundary and distance data
  curr_dist <- dist_list[[ct]]
  
  # Merge into metadata
  curr_obj$dist_to_boundary <- curr_dist$dist_to_boundary[
    match(row.names(curr_obj@meta.data), curr_dist$cellid)
  ]
  
  curr_obj$Cell_ID <- row.names(curr_obj@meta.data)
  
  all_meta_list[[ct]] <- curr_obj@meta.data
  assign(obj_name, curr_obj)
}

# Combine the meta data
final_meta_all <- bind_rows(all_meta_list, .id = "cancertype")


# (5). Merge seurat objects from all cancer types
ls()[grep("gradient_caners",ls())]
obj_names <- ls(pattern = "^obj_.*_gradient_caners$")
obj_names
obj_list <- mget(obj_names)
obj_all_gradient_cancers <- merge(x = obj_list[[1]], 
                                  y = obj_list[-1])

obj_all_gradient_cancers <- JoinLayers(obj_all_gradient_cancers)
obj_all_gradient_cancers

identical(row.names(obj_all_gradient_cancers@meta.data),row.names(final_meta_all))



################################################################################
# 1_2. Calculate gradient genes
################################################################################

# (1). Set parameters
bin_width <- 20
max_dist <- 200
filtered_meta <- final_meta_all %>%
  filter(sample %in% gradient_caners) %>%
  filter(!is.na(dist_to_boundary)) %>%
  filter(abs(dist_to_boundary) < 400) %>%
  mutate(dist_bin = floor(dist_to_boundary / bin_width) * bin_width)
unique(filtered_meta$sample)
head(filtered_meta)

exclude_types <- c("Keratinocyte","Oligodendrocyte")
target_cancertype <- unique(filtered_meta$cancertype)
target_celltypes <- setdiff(unique(filtered_meta$celltype), exclude_types)
target_cancertype
target_celltypes

# (2). Get data in 200um region
meta_200 <- filtered_meta %>%
  filter(abs(dist_to_boundary) <= max_dist) %>%
  mutate(bin_id = floor(dist_to_boundary / bin_width) * bin_width)

head(meta_200)
# Only data in tumor core (SN4)
meta_200 <- meta_200[!(meta_200$celltype == "Malignant/Epithelial" & meta_200$niche %in% c(1,2,3,5)),]

# Filter expr data in 200 um region
obj_200 <- subset(obj_all_gradient_cancers, cells = meta_200$Cell_ID)
exp_matrix <- GetAssayData(obj_200, assay = "Xenium", layer = "data")

# (3). Gradient gene calculation
directions <- list("Inward" = c(-200, 0), "Outward" = c(0, 200))
all_gradient_results_list <- list()

for (pri in target_cancertype) {
  message(">>> Analyzing Cancer Type: ", pri)
  for (ct in target_celltypes) {
    
    # Directions
    for (dir_name in names(directions)) {
      range <- directions[[dir_name]]
      
      # Specific tumor + specific cells + specific directional region
      current_meta <- meta_200 %>% 
        filter(cancertype == pri & 
                 celltype == ct & 
                 dist_to_boundary >= range[1] & 
                 dist_to_boundary <= range[2])
      
      # Cell count and Bin number
      if (nrow(current_meta) < 3) next
      unique_bins <- sort(unique(current_meta$bin_id))
      if (length(unique_bins) < 4) next # At least 4 points are needed to calculate the correlation.
      
      # Calculate the Bin mean matrix
      bin_mtx <- sapply(unique_bins, function(b) {
        cells <- current_meta$Cell_ID[current_meta$bin_id == b]
        if(length(cells) >= 3) {
          return(Matrix::rowMeans(exp_matrix[, cells, drop = FALSE]))
        } else {
          return(rep(NA, nrow(exp_matrix)))
        }
      })
      colnames(bin_mtx) <- unique_bins
      
      # Remove invalid lines
      bin_mtx <- bin_mtx[rowSums(!is.na(bin_mtx)) > (ncol(bin_mtx)/2), , drop=FALSE]
      if (nrow(bin_mtx) == 0) next
      
      # Calculate the correlation
      dist_vec <- as.numeric(colnames(bin_mtx))
      res <- apply(bin_mtx, 1, function(gene_row) {
        if(sd(gene_row, na.rm = TRUE) == 0) return(c(cor = 0, p = 1))
        
        test <- cor.test(dist_vec, gene_row, method = "pearson", use = "complete.obs")
        return(c(corr = as.numeric(test$estimate), p = as.numeric(test$p.value)))
      })
      
      # Summary of Results
      res_df <- as.data.frame(t(res))
      colnames(res_df) <- c("corr", "p")
      res_df <- res_df %>%
        mutate(gene = rownames(bin_mtx), 
               cancertype = pri, 
               celltype = ct,
               direction = dir_name,  # Direction
               adj_p = p.adjust(p, method = "fdr"))
      
      all_gradient_results_list[[paste(pri, ct, dir_name, sep = "_")]] <- res_df
    }
  }
}

# (4). Merge final result table and export
final_gradient_all <- bind_rows(all_gradient_results_list)
head(final_gradient_all)
write.csv(final_gradient_all, "~/space/gradient/1_1_Spatial_Gradient_bycancertype_all.csv")

final_gradient_sig <- final_gradient_all[final_gradient_all$p<0.05,]
dim(final_gradient_sig)
table(final_gradient_sig$cancertype)


################# 2. Quality control of candidate gradient genes
#################
#################
#################
#################
#################
#################


################################################################################
# 2_1. Calculate the gene expression proportion in each cell type across cancer types
################################################################################

library(dplyr)
library(Matrix)

cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types) {
  # Load data
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS files of each cancer were generated in Step02, line 135
  
  # Calculate the gene expression proportion in each cell type
  expr_matrix <- GetAssayData(obj, layer = "data")
  unique_groups <- unique(obj$celltype)
  all_groups_pct <- lapply(unique_groups, function(gp) {
    cells_in_gp <- colnames(obj)[obj$celltype == gp]
    
    if(length(cells_in_gp) == 0) return(NULL)
    pct_values <- rowSums(expr_matrix[, cells_in_gp, drop = FALSE] > 0) / length(cells_in_gp)
    group_info <- strsplit(gp, "-")[[1]]
    
    data.frame(
      gene = names(pct_values),
      pct = as.numeric(pct_values),
      celltype = group_info[1],
      cancertype = unique(obj$cancertype),
      stringsAsFactors = FALSE
    )
  })
  # Combine results
  df_all_pct_cancertype <- bind_rows(all_groups_pct)
  head(df_all_pct_cancertype)
  
  # Export to rds
  saveRDS(df_all_pct_cancertype,file=paste0("~/space/gradient/2_1_df_all_pct_celltype_",cancer,".rds"))
  
  
}

# Combine gene expression proportion data of every cancer types
df_all_pct <- dir("~/space/gradient/")[grep("2_1_df_all_pct_celltype_.*.rds",dir("~/space/gradient/"))]
df_all_pct <- list()
for (i in cancer_types) {
  df_all_pct[[i]] <- readRDS(paste0("~/space/gradient/2_1_df_all_pct_celltype_",i,".rds"))
}
df_all_pct <- do.call(rbind,df_all_pct)
head(df_all_pct)
saveRDS(df_all_pct, "~/space/gradient/2_1_df_all_pct_celltype_cancertype.rds")

################################################################################
# 2_2. Filter with gene expression proporition, perform filter by cancer type
################################################################################

# Perform filter
final_gradient_sig_propfilter <- final_gradient_sig %>%
  left_join(df_all_pct, by = c("gene", "celltype", "cancertype")) %>%
  filter(pct > 0.05) %>%
  dplyr::select(-pct)

# Check results and export
head(final_gradient_sig_propfilter)
dim(final_gradient_sig) # Raw counts
dim(final_gradient_sig_propfilter)            # Counts after filter
table(final_gradient_sig_propfilter$celltype)
write.csv(final_gradient_sig_propfilter, "~/space/gradient/2_2_gradient_sig_propfilter.csv")



################################################################################
# 2_3. Filter with gene expression values in tumor cells and in all the cell of the boundary region
################################################################################
# (1). Calculate the expression value of each gene at a granularity of cancer type × cell type × distance gradient (Bin).

# Set parameters
max_dist <- 800
meta_800 <- final_meta_all %>%
  filter(abs(dist_to_boundary) <= max_dist) %>%
  mutate(bin_id = floor(dist_to_boundary / bin_width) * bin_width)

# Tumor cells only in tumor core
meta_800 <- meta_800[!(meta_800$celltype == "Malignant/Epithelial" & meta_800$niche %in% c(1,2,3,5)),]
# Expr data
exp_matrix_800 <- GetAssayData(subset(obj_all_gradient_cancers, cells = row.names(meta_800)), assay = "Xenium", layer = "data")

# Calculte the expression value
directions <- list("Inward" = c(-800, 0), "Outward" = c(0, 800))

all_expression_long_list <- list()
for (samp in cancer_types) {
  message(">>> Processing Cancer Type: ", samp)
  for (ct in target_celltypes) {
    for (dir_name in names(directions)) {
      range <- directions[[dir_name]]
      
      # Metadata
      current_meta <- meta_800 %>% 
        filter(cancertype == samp & celltype == ct & 
                 dist_to_boundary >= range[1] & dist_to_boundary <= range[2])
      
      if (nrow(current_meta) < 3) next
      unique_bins <- sort(unique(current_meta$bin_id))
      if (length(unique_bins) < 4) next
      
      # Calculate the Bin mean matrix
      bin_mtx <- sapply(unique_bins, function(b) {
        cells <- current_meta$Cell_ID[current_meta$bin_id == b]
        if(length(cells) >= 3) {
          return(Matrix::rowMeans(exp_matrix_800[, cells, drop = FALSE]))
        } else {
          return(rep(NA, nrow(exp_matrix_800)))
        }
      })
      
      rownames(bin_mtx) <- rownames(exp_matrix_800)
      colnames(bin_mtx) <- unique_bins
      
      # Keep only rows that are not all NA
      valid_rows <- rowSums(!is.na(bin_mtx)) > 0
      
      if(any(valid_rows)){
        # Extract valid data
        sub_mtx <- bin_mtx[valid_rows, , drop = FALSE] 
        
        # Convert into long format
        temp_long <- as.data.frame(sub_mtx) %>%
          # Set gene name column
          tibble::rownames_to_column(var = "gene") %>%
          tidyr::pivot_longer(
            cols = -gene, 
            names_to = "bin_id", 
            values_to = "mean_expression"
          ) %>%
          mutate(
            cancertype= samp, 
            celltype = ct, 
            direction = dir_name,
            bin_id = as.numeric(bin_id)
          )
        
        all_expression_long_list[[paste(samp, ct, dir_name, sep = "_")]] <- temp_long
      }
    }
  }
}

# Merge into the final large expression table
final_gradient_expression_sample_df <- data.table::rbindlist(all_expression_long_list)
saveRDS(final_gradient_expression_sample_df, "~/space/gradient/2_3_gradient_expression_df_cancertype.rds")


# (2). Calculate the expression value of each gene at a granularity of sample  × distance gradient (Bin).

# Set parameters
directions <- list("Inward" = c(-800, 0), "Outward" = c(0, 800))
all_samples_allcells_list <- list()

# Calculation
for (samp in cancer_types) {
  message(">>> Processing All_Cells for Sample: ", samp)
  
  for (dir_name in names(directions)) {
    range <- directions[[dir_name]]
    
    current_meta <- meta_800 %>% 
      filter(cancertype == samp & dist_to_boundary >= range[1] & dist_to_boundary <= range[2])
    
    if (nrow(current_meta) < 3) next
    unique_bins <- sort(unique(current_meta$bin_id))
    if (length(unique_bins) < 4) next
    
    bin_mtx <- sapply(unique_bins, function(b) {
      cells <- current_meta$Cell_ID[current_meta$bin_id == b]
      if(length(cells) >= 3) {
        return(Matrix::rowMeans(exp_matrix_800[, cells, drop = FALSE]))
      } else {
        return(rep(NA, nrow(exp_matrix_800)))
      }
    })
    
    colnames(bin_mtx) <- unique_bins
    rownames(bin_mtx) <- rownames(exp_matrix_800)
    
    # Conver into long format
    temp_long <- as.data.frame(bin_mtx) %>%
      tibble::rownames_to_column(var = "gene") %>%
      tidyr::pivot_longer(
        cols = -gene, 
        names_to = "bin_id", 
        values_to = "mean_expression"
      ) %>%
      mutate(
        cancertype = samp, 
        celltype = "All_Cells", 
        direction = dir_name,
        bin_id = as.numeric(bin_id)
      )
    
    all_samples_allcells_list[[paste(samp, dir_name, sep = "_")]] <- temp_long
  }
}

# Merge and export
final_all_cells_expression_df <- data.table::rbindlist(all_samples_allcells_list)
head(final_all_cells_expression_df)
saveRDS(final_all_cells_expression_df, "~/space/gradient/2_3_gradient_expression_df_allcells.rds")


# (3). Filter genes whose max expression in tumor cells greater than max expression in the boundary region

# Calculate the maximum value of tumor cells in the range [-200, -20]
inward_max_df <- final_gradient_expression_sample_df %>%
  filter(bin_id >= -200 & bin_id <= -20 & celltype == "Malignant/Epithelial" & direction =="Inward") %>%
  group_by(cancertype, gene) %>%
  summarise(max_inward_expr = max(mean_expression, na.rm = TRUE), .groups = 'drop')
head(inward_max_df)

# Calculate the values in the boundary region
outward_val_df <- final_all_cells_expression_df %>%
  filter(bin_id %in% c(-40, -20, 0, 20, 40)) %>%
  group_by(cancertype, gene) %>%
  summarise(
    val_minus_40 = mean(mean_expression[bin_id == -40], na.rm = TRUE),
    val_minus_20 = mean(mean_expression[bin_id == -20], na.rm = TRUE),
    val_0        = mean(mean_expression[bin_id == 0], na.rm = TRUE),
    val_20       = mean(mean_expression[bin_id == 20], na.rm = TRUE),
    val_40       = mean(mean_expression[bin_id == 40], na.rm = TRUE),
    .groups = 'drop'
  )

# Merge and filter records that meet the criteria
result <- inward_max_df %>%
  inner_join(outward_val_df, by = c("cancertype", "gene")) %>%
  # The maximum value is greater than the maximum value in the boundary region (-20 - 40)
  filter(max_inward_expr > val_0 & max_inward_expr > val_20 & max_inward_expr > val_minus_20 )
head(result)

# (4). Perform filter
# Ensure that the cancertype and gene columns in the result match the column names in the list
filter_keys <- result %>% dplyr::select(cancertype, gene)

# Set the corr to -corr, as the distance inward tumor are negative
final_gradient_sig_propfilter$corr[final_gradient_sig_propfilter$direction=="Inward"]= -final_gradient_sig_propfilter$corr[final_gradient_sig_propfilter$direction=="Inward"] 
# Extract results of Epi and inward and negative
final_gradient_sig_propfilter_epi_inward_nega <- final_gradient_sig_propfilter[final_gradient_sig_propfilter$celltype=="Malignant/Epithelial" & final_gradient_sig_propfilter$direction=="Inward" & final_gradient_sig_propfilter$corr < 0,]
# Filter with max values
final_gradient_sig_propfilter_epi_inward_nega_filtered <-  final_gradient_sig_propfilter_epi_inward_nega %>% semi_join(filter_keys, by = c("cancertype", "gene"))
head(final_gradient_sig_propfilter_epi_inward_nega_filtered)


################################################################################
# 3_1. Cancer-shared gradient genes
################################################################################

# (1). Conserved negative gradient genes in tumor cells
head(final_gradient_sig_propfilter_epi_inward_nega_filtered)
final_gradient_sig_propfilter_epi_inward_nega_filtered$neg10p <- -log10(final_gradient_sig_propfilter_epi_inward_nega_filtered$p)

conserved_genes_table_neg_filter <- final_gradient_sig_propfilter_epi_inward_nega_filtered %>%
  filter(p < 0.05) %>%
  mutate(trend = ifelse(corr > 0, "Positive", "Negative")) %>%
  group_by(direction, celltype, gene, trend) %>%
  summarise(
    n_cancers = n_distinct(cancertype),
    cancer_list = paste(sort(unique(cancertype)), collapse = ", "),
    avg_corr = mean(corr),
    avg_log10p = mean(neg10p),
    .groups = "drop"
  ) %>%
  filter(n_cancers >= 1) %>%
  arrange(direction, celltype, desc(n_cancers), desc(abs(avg_corr)))

table(conserved_genes_table_neg_filter$celltype)
head(conserved_genes_table_neg_filter)

# Genes conserved in at least 6 cancer types
table(conserved_genes_table_neg_filter$celltype[conserved_genes_table_neg_filter$n_cancers>=6])
# Export
write.csv(conserved_genes_table_neg_filter, "~/space/gradient/3_1_conserved_genes_table_neg_filter.csv")


# (2). Conserved positive gradient genes in tumor cells
final_gradient_sig_propfilter$neg10p <- -log10(final_gradient_sig_propfilter$p)
conserved_genes_table_propfilter_all <- final_gradient_sig_propfilter %>%
  filter(p < 0.05) %>%
  mutate(trend = ifelse(corr > 0, "Positive", "Negative")) %>%
  group_by(direction, celltype, gene, trend) %>%
  summarise(
    n_cancers = n_distinct(cancertype),
    cancer_list = paste(sort(unique(cancertype)), collapse = ", "),
    avg_corr = mean(corr),
    avg_log10p = mean(neg10p),
    .groups = "drop"
  ) %>%
  filter(n_cancers >= 1) %>%
  arrange(direction, celltype, desc(n_cancers), desc(abs(avg_corr)))

head(conserved_genes_table_propfilter_all)

# Extract conserved positive gradient genes in tumor cells
conserved_genes_table_propfilter_epi_inward_pos <- conserved_genes_table_propfilter_all[conserved_genes_table_propfilter_all$direction=="Inward" & conserved_genes_table_propfilter_all$celltype=="Malignant/Epithelial" & conserved_genes_table_propfilter_all$trend=="Positive",]
head(conserved_genes_table_propfilter_epi_inward_pos)


# Combine negative and positive results
conserved_genes_table_filter <- rbind(conserved_genes_table_neg_filter, conserved_genes_table_propfilter_epi_inward_pos)
head(conserved_genes_table_filter)
tail(conserved_genes_table_filter)
write.csv(conserved_genes_table_filter, "~/space/gradient/3_1_conserved_genes_table_filter.csv")



################################################################################
# 4_1. Reactome enrichment based-on Cancer-shared gradient genes
################################################################################
library(clusterProfiler)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(biomaRt)
library(ReactomePA)

data_reactome_merge <- data.frame()
conserved_genes_table_filter_min6cancers <- conserved_genes_table_filter[conserved_genes_table_filter$n_cancers>5,]

# (1). Set parameters
celltype_list <- unique(conserved_genes_table_filter_min6cancers$celltype)
direction_list <- unique(conserved_genes_table_filter_min6cancers$direction)
trend_list <- unique(conserved_genes_table_filter_min6cancers$trend)

# (2). Perform Reactome enrichment
for (cell in celltype_list) {
  for (dir in direction_list) {
    for (tr in trend_list) {
      
      # Filter genes in current condition
      current_genes_df <- conserved_genes_table_ownpropfilter0.05_min6cancers %>%
        filter(celltype == cell & direction == dir & trend == tr)
      
      if (nrow(current_genes_df) == 0) next
      gene_symbols <- unique(current_genes_df$gene)
      
      tryCatch({
        # Gene id transform (Symbol -> ENTREZID)
        gene_ids <- bitr(gene_symbols, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
        
        # Enrichment
        reactome <- enrichPathway(gene = gene_ids$ENTREZID, 
                                  pvalueCutoff = 1, 
                                  qvalueCutoff = 1, 
                                  readable = TRUE)
        if (is.null(reactome) || nrow(as.data.frame(reactome)) == 0) {
          message(paste("No pathways enriched for:", cell, dir, tr))
          next
        }
        
        # Merge into a data.frame
        data_reactome <- as.data.frame(reactome)
        data_reactome$neglogpvalue <- -log10(data_reactome$pvalue)
        data_reactome$neglogqvalue <- -log10(data_reactome$qvalue)
        
        # Add labels
        data_reactome$celltype <- cell
        data_reactome$direction <- dir
        data_reactome$trend <- tr
        
        # Merge into the result table
        data_reactome_merge <- rbind(data_reactome_merge, data_reactome)
        
        message(paste("Processed Success:", cell, "-", dir, "-", tr, ":", nrow(data_reactome), "Pathways."))
        
      }, error = function(e) {
        message(paste("Error in", cell, dir, tr, ":", e$message))
      })
    }
  }
}

# Filter pathways with pvalue < 0.05
data_reactome_merge <- data_reactome_merge[data_reactome_merge$pvalue<0.05,]
data_reactome_merge <- data_reactome_merge[!is.na(data_reactome_merge$neglogqvalue),]

# (3). Export the results
write.csv(data_reactome_merge, "~/space/gradient/4_1_data_reactome_merge.csv", row.names = FALSE)





