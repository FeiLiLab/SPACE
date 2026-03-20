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
allcancer_meta_merged = allcancer_meta_merged[!(allcancer_meta_merged$celltype %in% c("Malignant/Epithelial",'Oligodendrocyte','Keratinocyte')),]
head(allcancer_meta_merged)

# (2). Split the metadata by spatial niche and calculated the cell proportion in each spatial niche
for (i in 1:5) {
  niche_df <- allcancer_meta_merged[allcancer_meta_merged$niche == i, ]
  assign(paste0("meta_niche", i), niche_df, envir = .GlobalEnv)
  
  tab <- table(niche_df$cancertype, niche_df$celltype) %>% as.data.frame.matrix()
  row_sums <- rowSums(tab)
  prop <- tab / row_sums
  prop_t <- t(prop)
  
  assign(paste0("meta_niche", i, "_celltype_percent"), as.data.frame(prop_t), envir = .GlobalEnv)
}

# (3). Fill NA with 0
setdiff(colnames(meta_niche1_celltype_percent),colnames(meta_niche5_celltype_percent))
meta_niche5_celltype_percent[,setdiff(colnames(meta_niche1_celltype_percent),colnames(meta_niche5_celltype_percent))]=0
meta_niche5_celltype_percent <- meta_niche5_celltype_percent[,colnames(meta_niche1_celltype_percent)]



# (4). Export the proportion matrix 
for (i in 1:5) {
  obj <- get(paste0("meta_niche",i,"_celltype_percent"))
  write.csv(obj,file=paste0("~/space/ORR/niche",i,"_celltype_percent.csv"))
}



################################################################################
# 1_2. Correlation analysis between cell type prop and objective response rate in cancer immune checkpoint therapy
################################################################################

# (1). Generate a matrix to store the results
correlation_of_cellprop_ORR=matrix(NA,nrow = 9,ncol=11)
row.names(correlation_of_cellprop_ORR)=rownames(meta_niche1_celltype_percent)
colnames(correlation_of_cellprop_ORR)=c("celltype",paste(rep(c("correlation_niche","p_value_niche"),5),rep(c(1:5),each=2),sep=""))
correlation_of_cellprop_ORR=as.data.frame(correlation_of_cellprop_ORR)
cancer_orr = data.frame("cancer"=c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD"),
                        "ORR" = c(23,7.1,4.3,15,8.7,15.2,17.03,0,5,19.9,37,20.4))
row.names(cancer_orr)=cancer_orr$cancer


# (2). Correlation analysis
for (j in 1:5) {
  obj <- get(paste0("meta_niche",j,"_celltype_percent"))
  for(i in 1:nrow(obj)) {
    
    row_data <- as.numeric(obj[i,])
    cancer_vector <- as.numeric(cancer_orr[colnames(obj), 2])
    cor_test <- cor.test(row_data, cancer_vector, method = "pearson")
    eval(parse(text = paste0("correlation_of_cellprop_ORR$correlation_niche",j,"[i] <- cor_test$estimate")))
    eval(parse(text = paste0("correlation_of_cellprop_ORR$p_value_niche",j,"[i] <- cor_test$p.value")))
    
  }
}
# Validation of correlation results in a single cell type
cor.test(as.numeric(meta_niche2_celltype_percent_mean[5,]),cancer_orr[, 2])

# (3). Correlation analysis (Remove the outlier SKCM) 
correlation_of_cellprop_ORR_noskcm <- correlation_of_cellprop_ORR
for (j in 1:5) {
  obj <- get(paste0("meta_niche",j,"_celltype_percent"))
  for(i in 1:nrow(obj)) {
    
    row_data <- as.numeric(obj[i,])
    cancer_vector <- as.numeric(cancer_orr[colnames(obj), 2])
    cor_test <- cor.test(row_data[-11], cancer_vector[-11], method = "pearson")
    eval(parse(text = paste0("correlation_of_cellprop_ORR_noskcm$correlation_niche",j,"[i] <- cor_test$estimate")))
    eval(parse(text = paste0("correlation_of_cellprop_ORR_noskcm$p_value_niche",j,"[i] <- cor_test$p.value")))
    
  }
}
cor.test(as.numeric(meta_niche2_celltype_percent_mean[5,-11]),cancer_orr[-11, 2])

# (4). Append cell type information to the results matrix
correlation_of_cellprop_ORR$celltype=row.names(correlation_of_cellprop_ORR)
correlation_of_cellprop_ORR_noskcm$celltype=row.names(correlation_of_cellprop_ORR_noskcm)

# (5). Export results
write.csv(correlation_of_cellprop_ORR,file="~/space/ORR/correlation_of_cellprop_ORR.csv",row.names = T)
write.csv(correlation_of_cellprop_ORR_noskcm,file="~/space/ORR/correlation_of_cellprop_ORR_noskcm.csv",row.names = T)

