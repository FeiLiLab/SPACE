library(dplyr)
library(data.table)
library(tidyverse)


################# 1. Correlation analysis between niche-resolved cell proportion and objective response rate (ORR) in cancer immune checkpoint therapy (ICB)
#################
#################
#################
#################
#################


################################################################################
# 1_1. Load metadata
################################################################################

allcancer_meta_merged = readRDS(file="~/space/seurat_niche/metadata_merge.rds") # Merged meta data of all cancers was generated in Step02, line 139

# (1). Remove cancer specify cell types and malignant
allcancer_meta_merged = allcancer_meta_merged[!(allcancer_meta_merged$celltype %in% c("Malignant",'Oligodendrocyte','Keratinocyte')),]
head(allcancer_meta_merged)

# (2). Split the metadata by spatial niche and calculated the cell proportion in each spatial niche
for (i in 1:5) {
  niche_df <- allcancer_meta_merged[allcancer_meta_merged$niche == i, ]
  assign(paste0("meta_niche", i), niche_df, envir = .GlobalEnv)
  
  tab <- table(niche_df$sample, niche_df$celltype) %>% as.data.frame.matrix()
  row_sums <- rowSums(tab)
  prop <- tab / row_sums
  prop_t <- t(prop)
  
  assign(paste0("meta_niche", i, "_celltype_percent"), prop_t, envir = .GlobalEnv)
}

# (3). Fill missing proportions with 0 when a sample has no cells in a given spatial niche
all_samples <- unique(allcancer_meta_merged$sample)
niches <- 1:5
for (i in niches) {
  mat_name <- paste0("meta_niche", i, "_celltype_percent")
  
  # Get data of current niche
  current_mat <- get(mat_name, envir = .GlobalEnv)
  
  # Get Samples in the current niche
  current_samples <- colnames(current_mat)
  
  # Get samples missing cells in the current niche
  missing_samples <- setdiff(all_samples, current_samples)
  
  # Fill cel proportion with 0
  if (length(missing_samples) > 0) {
    message(" Found ", length(missing_samples), " missing samples, fill with 0")
    
    # Create a matrix of 0
    zero_mat <- matrix(0,
                       nrow = nrow(current_mat),
                       ncol = length(missing_samples))
    
    colnames(zero_mat) <- missing_samples
    rownames(zero_mat) <- rownames(current_mat)
    
    # Combine the 0 matrix and current data
    current_mat <- cbind(current_mat, zero_mat)
    
    # Set the order of result matrix
    current_mat <- current_mat[, all_samples, drop = FALSE]
    
    # Assign the matrix to the raw name
    assign(mat_name, current_mat, envir = .GlobalEnv)
    
    message("Matrix filled, column number：", ncol(current_mat))
  } else {
    message("This niche don't need fill")
  }
}


# (4). Calculate the mean cell proportion in cancer level 
for (i in 1:5) {
  long_data <- as.data.frame(get(paste0("meta_niche",i,"_celltype_percent"))) %>%
    rownames_to_column("Feature") %>%
    pivot_longer(
      cols = -Feature,
      names_to = "Sample",
      values_to = "Value"
    )
  long_data$CancerType <- gsub("_.*", "",long_data$Sample)
  
  grouped_means <- long_data %>%
    group_by(Feature, CancerType) %>%
    summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  final_result <- grouped_means %>%
    pivot_wider(
      names_from = CancerType,
      values_from = Mean_Value
    ) %>%
    column_to_rownames("Feature") %>%
    as.data.frame()
  
  assign(paste0("meta_niche",i,"_celltype_percent_mean"),final_result)
}

# (5). Export the proportion matrix 
for (i in 1:5) {
  obj <- get(paste0("meta_niche",i,"_celltype_percent_mean"))
  write.csv(obj,file=paste0("~/space/ORR/niche",i,"_celltype_percent_mean_fill.csv"))
}



################################################################################
# 1_2. Correlation analysis between cell type prop and objective response rate in cancer immune checkpoint therapy
################################################################################

# (1). Generate a matrix to store the results
correlation_of_cellprop_ORR=matrix(NA,nrow = 9,ncol=11)
row.names(correlation_of_cellprop_ORR)=rownames(meta_niche1_celltype_percent_mean)
colnames(correlation_of_cellprop_ORR)=c("celltype",paste(rep(c("correlation_niche","p_value_niche"),5),rep(c(1:5),each=2),sep=""))
correlation_of_cellprop_ORR=as.data.frame(correlation_of_cellprop_ORR)
cancer_orr = data.frame("cancer"=c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD"),
                        "ORR" = c(23,7.1,4.3,15,8.7,15.2,17.03,0,5,19.9,37,20.4))
row.names(cancer_orr)=cancer_orr$cancer

# (2). Correlation test
for (j in 1:5) {
  obj <- get(paste0("meta_niche",j,"_celltype_percent_mean"))
  for(i in 1:nrow(obj)) {
    
    row_data <- as.numeric(obj[i,])
    cancer_vector <- as.numeric(cancer_orr[colnames(obj), 2])
    cor_test <- cor.test(row_data, cancer_vector, method = "pearson")
    eval(parse(text = paste0("correlation_of_cellprop_ORR$correlation_niche",j,"[i] <- cor_test$estimate")))
    eval(parse(text = paste0("correlation_of_cellprop_ORR$p_value_niche",j,"[i] <- cor_test$p.value")))
    
  }
}
# Validation of correlation results in a single cell type
cor.test(as.numeric(meta_niche1_celltype_percent_mean[1,]),cancer_orr[, 2])

# (3). Append cell type information to the results matrix
correlation_of_cellprop_ORR$celltype=row.names(correlation_of_cellprop_ORR)

# (4). Export results
dir.create("~/space/ORR/")
write.csv(correlation_of_cellprop_ORR,file="~/space/ORR/correlation_of_cellprop_ORR.csv",row.names = T)


