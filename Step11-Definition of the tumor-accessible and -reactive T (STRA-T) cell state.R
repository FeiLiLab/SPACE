library(Seurat)
library(dplyr)
library(BiocParallel)
library(SeuratWrappers)
library(future)
library(AUCell)
library(GSEABase)
library(msigdbr)
library(KEGGREST)
library(data.table)
options(future.globals.maxSize = 300 * 1024^3)



################# 1. Identify the solid tumor-accessible and -reactive T (STAR-T) cell state
#################
#################
#################
#################
#################
#################

# (1). Load public signature data of T cell state
# The signature data of T cell state was obtained from the publication: 'Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance'
tcellstateNatMed <- read.csv("~/space/analysis/signature/tcell/signature_tcellstate_NatMed2023_merge.csv") %>%
  split(.$Cellstate) %>%
  lapply(function(x) unique(x$Gene))
names(tcellstateNatMed)

# (2). Perform T cell state signature evaluation in all the cancer types
plan("multisession",workers=8)
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types) {
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds"))
  obj_select<-subset(obj,celltype=="T_cell" )
  # Extract expression matrix
  expr_matrix <- as.matrix(GetAssayData(obj_select, slot = "data")) 
  set.seed(123)
  # Run AUCell_buildRankings
  cells_rankings <- AUCell_buildRankings(
    expr_matrix, 
    plotStats = FALSE,  # ?ر?ͳ??ͼ????
    nCores = 8)
  # Run AUCell_calcAUC
  cells_AUC <- AUCell_calcAUC(
    tcellstateNatMed, 
    cells_rankings,
    nCores = 8)
  saveRDS(cells_AUC,paste0("~/space/analysis/domain/tcell/start/1_1_allTcells_AUC_tcellstateNatMed_",cancer,".rds"))
}

# (3). Construct signature score-cell matrix across all the cancer types
AUC_tcell_merge <- list()
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types) {
  cells_AUC <- readRDS(paste0("~/space/analysis/domain/tcell/start/1_1_allTcells_AUC_tcellstateNatMed_", cancer, ".rds"))
  cells_AUC_tcell <- as.data.frame(getAUC(cells_AUC))
  AUC_tcell_merge[[cancer]] <- cells_AUC_tcell
}
AUC_tcell_NatMed <- do.call(cbind, AUC_tcell_merge)
colnames(AUC_tcell_NatMed) <- gsub("^[a-z]+\\.", "", colnames(AUC_tcell_NatMed))

seurat_niches<-readRDS("~/space/seurat_niche/metadata_merge.rds") # Merged meta data of all cancers was generated in Step02, line 139
seurat_niches_tcell<-seurat_niches[which(seurat_niches$celltype=="T_cell"),]


dt_seurat <- as.data.table(seurat_niches_tcell, keep.rownames = "cell_id")
dt_auc <- as.data.table(as.matrix(t(AUC_tcell_NatMed)), keep.rownames = TRUE)
setnames(dt_auc, "rn", "cell_id")
df_AUC_tcell_NatMed <- merge(dt_seurat, dt_auc, by = "cell_id")
setDF(df_AUC_tcell_NatMed)
rownames(df_AUC_tcell_NatMed) <- df_AUC_tcell_NatMed$cell_id
df_AUC_tcell_NatMed$cell_id <- NULL

# (4). Save the result of signature scores
saveRDS(df_AUC_tcell_NatMed,"~/space/analysis/domain/tcell/start/1_2_df_AUC_tcell_NatMed.rds")
df_AUC_tcell_NatMed<-readRDS("~/space/analysis/domain/tcell/start/1_2_df_AUC_tcell_NatMed.rds")

# (5). Extract signature scores of cytotoxicity and exhaustion
niche_Cytotoxicity <- df_AUC_tcell_NatMed %>%
  group_by(niche) %>%
  summarise(
    mean_coord_x = mean(Cytotoxicity, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'  )
print(niche_Cytotoxicity)

niche_Exhaustion <- df_AUC_tcell_NatMed %>%
  group_by(niche) %>%
  summarise(
    mean_coord_x = mean(Exhaustion, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'  )
print(niche_Exhaustion)

# (6). Identify STAR-T cells
df_AUC_tcell_SN4<-df_AUC_tcell_NatMed[which(df_AUC_tcell_NatMed$niche=="4"),]

cytotoxicity_median <- median(df_AUC_tcell_SN4$Cytotoxicity[df_AUC_tcell_SN4$Cytotoxicity > 0], na.rm = TRUE)
exhaustion_median <- median(df_AUC_tcell_SN4$Exhaustion[df_AUC_tcell_SN4$Exhaustion > 0], na.rm = TRUE)

df_AUC_tcell_SN4_start <- df_AUC_tcell_SN4[
  df_AUC_tcell_SN4$Cytotoxicity > cytotoxicity_median & 
    df_AUC_tcell_SN4$Exhaustion < exhaustion_median,]

# (7). Save the results of STAR-T cells
cell_start<-row.names(df_AUC_tcell_SN4_start)
cell_other<-setdiff(row.names(df_AUC_tcell_SN4),cell_start)
df_AUC_tcell_SN4$group_start <- ifelse(row.names(df_AUC_tcell_SN4) %in% cell_start, "cell_start", "cell_other")

saveRDS(df_AUC_tcell_SN4,"~/space/analysis/domain/tcell/start/1_3_df_AUC_tcell_SN4.rds")
write.csv(df_AUC_tcell_SN4,"~/space/analysis/domain/tcell/start/1_3_df_AUC_tcell_SN4.csv")

