library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(BiocParallel)
library(arrow)
library(future)
options(future.globals.maxSize = 300 * 1024^3) 

################# 1. Perform AUCell to evaluate metabolic activity based on REACTOME dataset
#################
#################
#################
#################
#################
library(AUCell)
library(GSEABase)
library(msigdbr)
library(KEGGREST)


################################################################################
# 1_1. Load metabolic signature datasets
################################################################################

meta_reactome <- GSEABase::getGmt("~/space/analysis/signature/metabolism/REACTOME_metabolism.gmt") # REACTOME pathways were obtained from The Molecular Signatures Database (MSigDB)

################################################################################
# 1_2. Run AUCell to evaluate metabolic activity
################################################################################

plan("multisession",workers=8)
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types) {
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS files of each cancer were generated in Step02, line 135
  # Extract expression matrix
  expr_matrix <- as.matrix(GetAssayData(obj, slot = "data")) 
  set.seed(123)
  # Run AUCell_buildRankings
  cells_rankings <- AUCell_buildRankings(
    expr_matrix, 
    plotStats = FALSE,
    nCores = 8)
  # Run AUCell based on Reactome dataset
  cells_AUC_reactome <- AUCell_calcAUC(
    meta_reactome, 
    cells_rankings,
    nCores = 8)
  saveRDS(cells_AUC_reactome,paste0("~/space/analysis/domain/metabolism/1_2_cells_AUC_reactome_",cancer,".rds"))
}


################# 2. Identify significantly differentially activated metabolic pathways among spatial niches
#################
#################
#################
#################
#################


################################################################################
# 2_1. Load metadata information
################################################################################

seurat_niches<-readRDS("~/space/seurat_niche/metadata_merge.rds") # Merged meta data of all cancers was generated in Step02, line 139
head(seurat_niches)
table(seurat_niches$majorcelltype)

################################################################################
# 2_2. Identify significantly differentially activated metabolic pathways based on REACTOME dataset
################################################################################

cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
auc_reactome<-data.frame()
for (cancer in cancer_types) {
  seurat_niches_select <- seurat_niches[which(seurat_niches$cancertype == cancer), ]
  majorcelltype_list <- names(table(seurat_niches_select$majorcelltype))
  cells_AUC_reactome_select <- getAUC(readRDS(paste0("~/space/analysis/domain/metabolism/1_2_cells_AUC_reactome_",cancer,".rds")))
  for (majorcelltype in majorcelltype_list) {
    niche_list <- names(table(seurat_niches_select[which(seurat_niches_select$majorcelltype == majorcelltype), ]$niche))
    for (niche in niche_list) {
      cell_target <- row.names(seurat_niches_select[which(seurat_niches_select$majorcelltype == majorcelltype & seurat_niches_select$niche == niche), ])
      cell_other <- setdiff(row.names(seurat_niches_select[which(seurat_niches_select$majorcelltype == majorcelltype), ]), cell_target)
      # Add empty value check
      if (length(cell_target) == 0 || length(cell_other) == 0) {
        next  # Skip this case
      }
      data_target <- cells_AUC_reactome_select[, cell_target, drop = FALSE]
      data_other <- cells_AUC_reactome_select[, cell_other, drop = FALSE]
      # Check number of data rows
      if (nrow(data_target) == 0 || nrow(data_other) == 0) {
        next
      }
      # Define function to perform significant analysis
      calculate_stats <- function(i) {
        target_row <- as.numeric(data_target[i, ])
        other_row <- as.numeric(data_other[i, ])
        mean_target <- mean(target_row, na.rm = TRUE)
        mean_other <- mean(other_row, na.rm = TRUE)
        se_target <- sd(target_row, na.rm = TRUE) / sqrt(sum(!is.na(target_row)))
        se_other <- sd(other_row, na.rm = TRUE) / sqrt(sum(!is.na(other_row)))
        if (mean_other == 0) {
          log2fc <- ifelse(mean_target == 0, 0, 20 * sign(mean_target))
        } else {
          log2fc <- log2((mean_target + 1e-10) / (mean_other + 1e-10))
        }
        
        p_val <- tryCatch({
          wilcox.test(target_row, other_row, exact = FALSE, na.rm = TRUE)$p.value
        }, error = function(e) NA)
        
        return(c(mean_target,mean_other,se_target,se_other,log2fc,p_val))
      }
      # Perform significant analysis
      # Use safe number of rows
      n_rows <- nrow(data_target)
      if (n_rows > 0) {
        stats_list <- lapply(1:n_rows, calculate_stats)
        results_reactome <- do.call(rbind, stats_list)
        
        df_reactome <- data.frame(
          pathway_name = row.names(data_target),
          cancertype = rep(cancer, n_rows),
          majorcelltype = rep(majorcelltype, n_rows),
          niche = rep(niche, n_rows),
          mean_target = results_reactome[, 1],
          mean_other = results_reactome[, 2],
          se_target = results_reactome[, 3],
          se_other = results_reactome[, 4],
          log2fc = results_reactome[, 5],
          p_value = results_reactome[, 6]
        )
        df_reactome$pathway_name <- gsub('["\']', '', df_reactome$pathway_name)
        df_reactome$p_adj <- p.adjust(df_reactome$p_value, method = "BH")
        df_reactome$change <- ifelse(df_reactome$p_value < 0.05 & df_reactome$log2fc > 0, "up",
                                     ifelse(df_reactome$p_value < 0.05 & df_reactome$log2fc < 0, "down", "nonsignificant"))
        auc_reactome<-rbind(auc_reactome,df_reactome)
        message("Processed: ", cancer,"-",majorcelltype,"-Niche",niche)
      }
    }
  }
}
write.csv(auc_reactome,"~/space/analysis/domain/metabolism/2_2_auc_reactome.csv")

################# 3. Identify cross-cancer shared dysregulated metabolic pathways
#################
#################
#################
#################
#################


################################################################################
# 3_1. Identify cross-cancer shared up-regulated metabolic pathways based on REACTOME dataset
################################################################################

auc_reactome<-read.csv("~/space/analysis/domain/metabolism/2_3_auc_reactome.csv")

auc_reactome_up<-auc_reactome[which(auc_reactome$change=="up"),]
majorcelltype_list <- names(table(auc_reactome_up$majorcelltype))
crosscancer_reactome_up<-data.frame()
for (majorcelltype in majorcelltype_list) {

  auc_reactome_select<-auc_reactome_up[which(auc_reactome_up$majorcelltype==majorcelltype),]
  niche_list<-names(table(auc_reactome_select$niche))
  for (niche in niche_list){

    pathway_means <- as.data.frame(auc_reactome_select %>%
                                     group_by(pathway_name) %>%
                                     summarise(
                                       mean_log2fc = mean(log2fc, na.rm = TRUE),
                                       mean_mean_target = mean(mean_target, na.rm = TRUE)))
    row.names(pathway_means)<-pathway_means$pathway_name
    df_crosscancer_reactome<-data.frame(table(auc_reactome_select[which(auc_reactome_select$niche==niche),]$pathway_name))
    df_crosscancer_reactome$majorcelltype<-rep(majorcelltype,nrow(df_crosscancer_reactome))
    df_crosscancer_reactome$niche<-rep(niche,nrow(df_crosscancer_reactome))
    df_crosscancer_reactome$mean_log2fc<-pathway_means[df_crosscancer_reactome$Var1,]$mean_log2fc
    df_crosscancer_reactome$mean_mean_target<-pathway_means[df_crosscancer_reactome$Var1,]$mean_mean_target
    crosscancer_reactome_up<-rbind(crosscancer_reactome_up,df_crosscancer_reactome)
    message("Processed: ", majorcelltype,"-Niche",niche)
  }
}
colnames(crosscancer_reactome_up)<-c("pathway_name","frequency","majorcelltype","niche","mean_log2fc","mean_mean_target")
write.csv(crosscancer_reactome_up,"~/space/analysis/domain/metabolism/3_1_crosscancer_reactome_up.csv")


################################################################################
# 3_2. Identify cross-cancer shared up-regulated metabolic pathways based on REACTOME dataset
################################################################################

auc_reactome_down<-auc_reactome[which(auc_reactome$change=="down"),]
majorcelltype_list <- names(table(auc_reactome_down$majorcelltype))
crosscancer_reactome_down<-data.frame()
for (majorcelltype in majorcelltype_list) {

  auc_reactome_select<-auc_reactome_down[which(auc_reactome_down$majorcelltype==majorcelltype),]
  niche_list<-names(table(auc_reactome_select$niche))
  for (niche in niche_list){

    pathway_means <- as.data.frame(auc_reactome_select %>%
                                     group_by(pathway_name) %>%
                                     summarise(
                                       mean_log2fc = mean(log2fc, na.rm = TRUE),
                                       mean_mean_target = mean(mean_target, na.rm = TRUE)))
    df_crosscancer_reactome<-data.frame(table(auc_reactome_select[which(auc_reactome_select$niche==niche),]$pathway_name))
    df_crosscancer_reactome$majorcelltype<-rep(majorcelltype,nrow(df_crosscancer_reactome))
    df_crosscancer_reactome$niche<-rep(niche,nrow(df_crosscancer_reactome))
    df_crosscancer_reactome$mean_log2fc<-pathway_means[df_crosscancer_reactome$Var1,]$mean_log2fc
    df_crosscancer_reactome$mean_mean_target<-pathway_means[df_crosscancer_reactome$Var1,]$mean_mean_target
    crosscancer_reactome_down<-rbind(crosscancer_reactome_down,df_crosscancer_reactome)
    message("Processed: ", majorcelltype,"-Niche",niche)
  }
}
colnames(crosscancer_reactome_down)<-c("pathway_name","frequency","majorcelltype","niche","mean_log2fc","mean_mean_target")
write.csv(crosscancer_reactome_down,"~/space/analysis/domain/metabolism/3_2_crosscancer_reactome_down.csv")


################# 4. Perform paired comparison to identify dysregulated metabolic pathways
#################
#################
#################
#################
#################


################################################################################
# 4_1. Perform paired comparison to identify dysregulated metabolic pathways based on REACTOME dataset
################################################################################

majorcelltype_list <- names(table(auc_reactome$majorcelltype))
paired_crosscancer_reactome<-data.frame()
for (majorcelltype in majorcelltype_list){
  auc_reactome_select1<-auc_reactome[which(auc_reactome$majorcelltype==majorcelltype),]
  niche_list<-names(table(auc_reactome_select1$niche))
  for (niche in niche_list){
    auc_reactome_select2<-auc_reactome_select1[which(auc_reactome_select1$niche==niche),]
    pathway_list<-names(table(auc_reactome_select2$pathway_name))
    for (pathway in pathway_list){
      auc_reactome_select3<-auc_reactome_select2[which(auc_reactome_select2$pathway_name==pathway),]
      p_value_select<-wilcox.test(auc_reactome_select3$mean_target,
                                  auc_reactome_select3$mean_other,
                                  paired = TRUE)[["p.value"]]
      log2fc_select<-log(mean(auc_reactome_select3$mean_target)/mean(auc_reactome_select3$mean_other),2)
      diff_select<-mean(auc_reactome_select3$mean_target)-mean(auc_reactome_select3$mean_other)
      crosscancer_reactome<-data.frame(pathway_name=pathway,
                                       majorcelltype=majorcelltype,
                                       niche=niche,
                                       p_value=p_value_select,
                                       log2fc=log2fc_select,
                                       n_cancer=nrow(auc_reactome_select3))
      paired_crosscancer_reactome<-rbind(paired_crosscancer_reactome,crosscancer_reactome)
    }
  }
}
paired_crosscancer_reactome$change <- ifelse(paired_crosscancer_reactome$p_value < 0.05 & paired_crosscancer_reactome$log2fc > 0, "up",
                                             ifelse(paired_crosscancer_reactome$p_value < 0.05 & paired_crosscancer_reactome$log2fc < 0, "down", "nonsignificant"))

write.csv(paired_crosscancer_reactome,"~/space/analysis/domain/metabolism/4_1_paired_crosscancer_reactome.csv")


