library(Seurat)
library(Matrix)
library(tidyr)
library(stringr)
library(dplyr)
library(tibble)
library(data.table)
library(tidyverse)
options(future.globals.maxSize = 300 * 1024^3)


################# 1. Identify differentially expressed genes (DEGs) between STAR-T and Other-T cells in spatial niche 4
#################
#################
#################
#################
#################
#################


################################################################################
# 1_1. Identify the DEGs between STAR-T and Other-T cells in spatial niche 4 regardless of cancer type
################################################################################

# (1). Load data
cancer_types <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (cancer in cancer_types) {
  obj <- readRDS(paste0("~/sapce/rds/obj_", cancer, ".rds"))
  assign(paste0("obj_",cancer), obj)
}

# (2). Merge data of T cells from all cancers into a combined Seurat object
tcells_list=list()
for (cancer in cancer_types) {
  tcells_list[[cancer]]=subset(get(paste0("obj_",cancer)), celltype== "T_cell")
}
tcell=merge(tcells_list[[1]],tcells_list[-1])
saveRDS(tcell,file="~/space/rds/tcell.rds")


# (3). Add cell subtype information to the Seurat object
df_AUC_tcell_SN4 <- readRDS("~/space/analysis/domain/tcell/start/1_3_df_AUC_tcell_SN4.rds") # Group information of STAR-T and other T cells was generated in Step11, line 111

cell_start <- row.names(df_AUC_tcell_SN4[df_AUC_tcell_SN4$group_start == "cell_start",])
cell_other <- row.names(df_AUC_tcell_SN4[df_AUC_tcell_SN4$group_start == "cell_other",])

tcell$group1 <- ifelse(colnames(tcell) %in% cell_start, "cell_start",
                       ifelse(colnames(tcell) %in% cell_other, "cell_other", "cell_otherniche"))

# (4). Perform DEG analysis between STAR-T vs. Other-T in spatial niche 4
tcell_subset1 <- subset(tcell, subset = group1 %in% c("cell_start", "cell_other"))
tcell_subset1@images=list()
tcell_subset1 <- JoinLayers(tcell_subset1)
Idents(tcell_subset1)<-tcell_subset1$group1
table(Idents(tcell_subset1))
deg_target<-FindMarkers(tcell_subset1, ident.1 = "cell_start",ident.2 = "cell_other")
deg_target$gene<-row.names(deg_target)

# (5). Perform DEG analysis between STAR-T vs. T cells in other spatial niches
tcell_subset2 <- subset(tcell, subset = group1 %in% c("cell_start", "cell_otherniche"))
tcell_subset2 <- JoinLayers(tcell_subset2)
Idents(tcell_subset2)<-tcell_subset2$group1
deg_other<-FindMarkers(tcell_subset2, ident.1 = "cell_start",ident.2 = "cell_otherniche")
deg_other$gene<-row.names(deg_other)

# (6). Save DEGs
write.csv(deg_target,"~/space/analysis/domain/tcell/start/1_4_deg_target.csv")
write.csv(deg_other,"~/space/analysis/domain/tcell/start/1_4_deg_other.csv")

deg_other_sig <- deg_other[deg_other$p_val_adj<0.05,]
deg_other_sig_up <- deg_other_sig[deg_other_sig$avg_log2FC>0,]
write.csv(deg_other_sig,"~/space/analysis/domain/tcell/start/1_4_deg_other_sig.csv")


# (7). Retain candidate genes that were consistently detected (Percentage > 5%) in the corresponding cell type in the public pan-cancer scRNA-seq dataset
# Gene expression percentage data were derived from the public pan-cancer scRNA-seq dataset reported in "Cross-tissue multicellular coordination and its rewiring in cancer"

pancancer_gene<-read.csv("~/space/analysis/signature/pancancer_celllineage.genes.percent.csv")

deg_target_sig <- deg_target[deg_target$p_val_adj < 0.05,]
deg_target_sig$neglogpadj <- ifelse(
  deg_target_sig$p_val_adj == 0,
  300,
  -log10(deg_target_sig$p_val_adj)
)

pancancer_gene_select <- pancancer_gene[which(pancancer_gene$Celltype == "T_cell" & pancancer_gene$Proportion > 0.05), ]
gene_expressed <- pancancer_gene_select$Gene 
deg_target_sig <- deg_target_sig[deg_target_sig$gene %in% gene_expressed, ]

deg_target_up <- deg_target_sig[deg_target_sig$avg_log2FC>0 & deg_target_sig$p_val_adj < 0.05,]
deg_target_down <- deg_target_sig[deg_target_sig$avg_log2FC<0 & deg_target_sig$p_val_adj < 0.05,]

write.csv(deg_target_sig,"~/space/analysis/domain/tcell/1_4_deg_target_sig.csv")


################################################################################
# 1_2. Identify cancer-shared DEGs in STAR-T cells
################################################################################

# (1). Identify of cancer-shared DEGs
deg_cancer_target <- data.frame()
for (cancer in unique(tcell_subset1$cancertype)) {
  obj <- subset(tcell_subset1, cancertype == cancer)
  deg_cancer_target_1 <- FindMarkers(obj, ident.1 = "cell_start",ident.2 = "cell_other")
  deg_cancer_target_1$cancer <- rep(cancer, nrow(deg_cancer_target_1))
  deg_cancer_target_1$gene <- row.names(deg_cancer_target_1)
  deg_cancer_target <- rbind(deg_cancer_target, deg_cancer_target_1)
}

# (2). Extract up-regulated genes in STAR-T cells
deg_cancer_target_up <- deg_cancer_target[deg_cancer_target$p_val_adj<0.05 & deg_cancer_target$avg_log2FC>0,]

# (3). Retain candidate genes that were consistently detected (Percentage > 5%) in the corresponding cell type in the public pan-cancer scRNA-seq dataset
deg_cancer_target_up <- deg_cancer_target_up[deg_cancer_target_up$gene %in% gene_expressed, ]

# (4). Retain candidate genes that were: (1) up-regulated in STAR-T versus other T cells in spatial niche 4 (pan-cancer), and (2) up-regulated in STAR-T versus T cells in the other four spatial niches (pan-cancer)
deg_cancer_target_up <- deg_cancer_target_up[which(deg_cancer_target_up$gene %in% deg_target_up$gene & deg_cancer_target_up$gene %in% deg_other_sig$gene[deg_other_sig$avg_log2FC>0]),]


# (5). Save results of cancer-shared DEGs in STAR-T cells
deg_cancer_target_up_final=data.frame(table(deg_cancer_target_up_final$gene)[order(table(deg_cancer_target_up_final$gene),decreasing = T)])
deg_cancer_target_up_final <- deg_cancer_target_up_final[deg_cancer_target_up_final$Freq>2,]
write.csv(deg_cancer_target_up,"~/space/analysis/domain/tcell/start/1_5_deg_cancer_target_up.csv")
write.csv(deg_cancer_target_up_final,"~/space/analysis/domain/tcell/start/1_5_df_deg_cancer_target_up_final.csv")



################# 2. Perform pathway enrichment analysis based on DEGs in STAR-T cells
#################
#################
#################
#################
#################
#################

library(clusterProfiler)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(biomaRt)
library(ReactomePA)

# (1). Perform pathway enrichment on up-regulated DEGs in STAR-T cells
geneid_up<- bitr(unique(deg_target_up$gene), fromType = "SYMBOL",toType = c("ENSEMBL", "SYMBOL","ENTREZID"),OrgDb = org.Hs.eg.db)
reactome_up <- enrichPathway(gene=geneid_up$ENTREZID,pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable=TRUE)
data.frame(reactome_up)$Description
df_reactome_up<-data.frame(reactome_up)
df_reactome_up$neglogqvalue<-(-log(df_reactome_up$qvalue,10))
df_reactome_up$Description<-factor(df_reactome_up$Description,levels=rev(df_reactome_up$Description))
write.csv(df_reactome_up,"~/space/analysis/domain/tcell/start/2_1_reactome_up.csv")

# (2). Perform pathway enrichment on down-regulated DEGs in STAR-T cells
geneid_down<- bitr(unique(deg_target_down$gene), fromType = "SYMBOL",toType = c("ENSEMBL", "SYMBOL","ENTREZID"),OrgDb = org.Hs.eg.db)
reactome_down <- enrichPathway(gene=geneid_down$ENTREZID,pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable=TRUE)
data.frame(reactome_down)$Description
df_reactome_down<-data.frame(reactome_down)
df_reactome_down$neglogqvalue<-(-log(df_reactome_down$qvalue,10))
df_reactome_down$Description<-factor(df_reactome_down$Description,levels=rev(df_reactome_down$Description))
write.csv(df_reactome_down,"~/space/analysis/domain/tcell/start/2_1_reactome_down.csv")


################# 3. Identify differential activated TFs of STAR-T compared with other T cells
#################
#################
#################
#################
#################
#################

################################################################################
# 3_1. Identify differential activated TFs (DATFs) in STAR-T vs. other T cells in SN4
################################################################################

# (1). Load results of pySCENIC TF analysis
TF_seurat_list <- saveRDS(file = "~/space/pyscenic/analysis/TF_RegulonAUC.rds") # RDS of pySCENIC TF analysis results was generated in Step06_1, line 194

# (2). Subset the Seurat object for T cells
tcell_list <- lapply(TF_seurat_list, `[[`, "T_cell")
seurat_tcell_merged <- merge(tcell_list[[1]], tcell_list[-1])
seurat_tcell_merged_meta.data = seurat_tcell_merged@meta.data

# (3). Load group information of STAR-T and other T cells
start_cells <- readRDS("~/space/analysis/domain/tcell/23_4_df_AUC_tcell_SN4.rds") # Group information of STAR-T and other T cells was generated in Step11, line 111
start_cells$cell_id <- as.character(start_cells$cell_id)
start_cells$group <- as.character(start_cells$group)

# (4). Map cell ids of STAR-T and other T cells
seurat_tcell_merged@meta.data$celltype[row.names(seurat_tcell_merged_meta.data) %in% start_cells$cell_id[start_cells$group=="STAR-T"]]="STAR-T"
seurat_tcell_merged@meta.data$celltype[row.names(seurat_tcell_merged_meta.data) %in% start_cells$cell_id[start_cells$group=="Other-T"]]="Other-T"
seurat_tcell_merged@meta.data$celltype[seurat_tcell_merged@meta.data$celltype == "T_cell"]="Tcell_otherniche"

unique(seurat_tcell_merged@meta.data$celltype)

# (5). Perform differential activated TFs analysis: STAR-T vs. Other-T in SN4
tcell_subset1 <- subset(seurat_tcell_merged, subset = celltype %in% c("STAR-T", "Other-T"))
Idents(tcell_subset1)<-tcell_subset1$celltype
datf_target<-FindMarkers(tcell_subset1, ident.1 = "STAR-T",ident.2 = "Other-T")
datf_target$gene<-row.names(datf_target)

# (6). Retain DATFs that were consistently detected (Percentage > 5%) in the corresponding cell type in the public pan-cancer scRNA-seq dataset
pancancer_gene_select<-pancancer_gene[which(pancancer_gene$majorcelltype=="T_cell"&pancancer_gene$Proportion>0.05),]
gene_expressed<-pancancer_gene_select$Gene 
datf_target<-datf_target[which(datf_target$gene %in% gene_expressed),]
datf_target$neglogpadj = -log10(datf_target$p_val_adj)
datf_target$gene = gsub("\\(\\+\\)","",datf_target$gene)
datf_target_sig=datf_target[datf_target$p_val_adj<0.05,]

datf_target_sig <- datf_target_sig %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))

# (7). Save the final results
write.csv(datf_target,"~/space/pyscenic/analysis/start/3_1_datf_start.csv")
write.csv(datf_target_sig,"~/space/pyscenic/analysis/start/3_1_datf_start_sig_ref0.05.csv")


################################################################################
# 3_2. Identify cancer-shared DATFs of STAR-T
################################################################################

# (1). Get all unique cancer types from metadata
cancer_types <- unique(tcell_subset1@meta.data$cancer)

# (2). Initialize a list to store DATF results for each cancer type
datf_list <- list()

# (3). Loop through each cancer type
for (cancer_type in cancer_types) {
  
  cat("Processing cancer:", cancer_type, "\n")
  
  # Subset T cells for the current cancer type
  tcell_current <- subset(tcell_subset1, subset = cancer == cancer_type)
  
  # Check if both STAR-T and Other-T groups exist to avoid FindMarkers errors
  if (length(unique(Idents(tcell_current))) < 2) {
    message("Skipping ", cancer_type, ": only one identity group present")
    next
  }
  
  # Calculate DATFs within the specific cancer type
  datf <- FindMarkers(tcell_current, 
                      ident.1 = "STAR-T", 
                      ident.2 = "Other-T")
  
  # If results are found, add gene names and cancer type labels
  if (nrow(datf) > 0) {
    datf$gene <- rownames(datf)
    datf$cancer <- cancer_type
    datf_list[[cancer_type]] <- datf
  }
}

# (4). Merge DATF results from all cancer types into a single matrix
datf_target_per_cancer <- bind_rows(datf_list, .id = "cancer")

# (5). Filter for significant DATFs
datf_target_per_cancer_sig <- datf_target_per_cancer[datf_target_per_cancer$p_val_adj < 0.05, ]

# (6). Calculate -log10(adjusted p-value) for downstream visualization
datf_target_per_cancer_sig$neglogpadj = -log10(datf_target_per_cancer_sig$p_val_adj)

# (7). Clean gene names
datf_target_per_cancer_sig$gene = gsub("\\(\\+\\)", "", datf_target_per_cancer_sig$gene)

# (8). Filter genes based on the public pan-cancer scRNA-seq dataset
pancancer_gene_select <- pancancer_gene[which(pancancer_gene$majorcelltype == "T_cell" & pancancer_gene$Proportion > 0.05), ]
gene_expressed <- pancancer_gene_select$Gene 
datf_target_per_cancer_sig <- datf_target_per_cancer_sig[which(datf_target_per_cancer_sig$gene %in% gene_expressed), ]

# (8). Export the significant DATF results
datf_target_per_cancer_up = datf_target_per_cancer_sig[datf_target_per_cancer_sig$avg_log2FC > 0, ]
datf_target_per_cancer_down = datf_target_per_cancer_sig[datf_target_per_cancer_sig$avg_log2FC < 0, ]
write.csv(datf_target_per_cancer_sig, "~/space/pyscenic/analysis/start/3_2_datf_start_per_cancer_sig_ref0.05.csv")

# (9). Calculate the frequency of cancers in which the TF is conserved and export results
TF_freq_up <- datf_target_per_cancer_up %>%
  group_by(gene) %>%
  summarise(Freq = n_distinct(cancer), .groups = "drop") %>%
  arrange(desc(Freq))
TF_freq_up=TF_freq_up[TF_freq_up$gene %in% datf_target_sig[datf_target_sig$direction == "Up", ]$gene,]
write.csv(TF_freq_up, "~/space/pyscenic/analysis/start/5_2_TF_frequency_across_cancers_up.csv", row.names = FALSE)

TF_freq_down <- datf_target_per_cancer_down %>%
  group_by(gene) %>%
  summarise(Freq = n_distinct(cancer), .groups = "drop") %>%
  arrange(desc(Freq))
TF_freq_down=TF_freq_down[TF_freq_down$gene %in% datf_target_sig[datf_target_sig$direction == "Down", ]$gene,]
write.csv(TF_freq_down, "~/space/pyscenic/analysis/start/3_2_TF_frequency_across_cancers_down.csv", row.names = FALSE)


################################################################################
# 3_3. Prepare data for visualization of cancer-shared DATFs in STAR-T
################################################################################

# (1). Remove suffixes (+) from the regulon name
extract_tf <- function(regulon_name) {
  gsub("\\(.+\\)", "", regulon_name) %>% 
    trimws() 
}

# (2). Collect all unique TF names across all tumors and cell types
all_tf_names <- character(0)

for (tumor in names(regulons_list)) {
  for (celltype in names(regulons_list[[tumor]])) {
    # Extract clean TF names for the current tumor/cell type subset
    tf_names <- extract_tf(names(regulons_list[[tumor]][[celltype]]))
    # Update the global list of unique TFs
    all_tf_names <- unique(c(all_tf_names, tf_names))
  }
}

# (3). Merge target genes for each TF across all tumors and cell types
merged_regulons <- list()

for (tf in all_tf_names) {
  
  # Temporary vector to accumulate target genes for the current TF
  all_targets <- character(0)
  
  for (tumor in names(regulons_list)) {
    for (celltype in names(regulons_list[[tumor]])) {
      
      # Get original regulon names for the current subset
      regulon_names <- names(regulons_list[[tumor]][[celltype]])
      
      # Identify regulons matching the current TF
      matching_idx <- which(extract_tf(regulon_names) == tf)
      
      if (length(matching_idx) > 0) {
        # Extract and aggregate target genes from all matching regulons
        targets_here <- unlist(regulons_list[[tumor]][[celltype]][matching_idx])
        all_targets <- c(all_targets, targets_here)
      }
    }
  }
  
  # Remove duplicates and save to the final list
  if (length(all_targets) > 0) {
    merged_regulons[[tf]] <- unique(all_targets)
  }
}
# (4). Save the final merged list to an RDS file
saveRDS(merged_regulons, "~/space/pyscenic/analysis/start/3_3_merged_regulons_cross_cancer.rds")



################# 4. Identify dysregulated ligand-receptor pairs in STAR-T cells
#################
#################
#################
#################
#################
#################

################################################################################
# 4_1. Calculate the cell-cell communication scores of STAR-T and Other-T
################################################################################

# (1). Load data
allcancer_niche4=readRDS("allcancer_niche4.rds") # RDS of cells in spatial niches was generated in Step08, line 62
allcancer_niche4 <- JoinLayers(allcancer_niche4)
start_cells <- readRDS("~/space/analysis/domain/tcell/23_4_df_AUC_tcell_SN4.rds") # Group information of STAR-T and other T cells was generated in Step11, line 111
CCC_data <- read.csv("~/space/ccc/cellchatdb_cpdb_italk.csv") 
# Human ligand-receptor interactions from CellChat (version 2.2.0), CellPhoneDB (version 5.0.0), and iTALK were merged, and only one-to-one L-R pairs were retained.

start_cells$cell_id <- as.character(start_cells$Cell_ID)
start_cells$group <- as.character(start_cells$group)

# (2). Map STAR-T cell id in metadata and Seurat object, set the group of T cell state
mapping_vector <- start_cells$group
names(mapping_vector) <- start_cells$cell_id
meta_cells <- rownames(allcancer_niche4@meta.data)
matched <- meta_cells %in% start_cells$cell_id
allcancer_niche4@meta.data$majorcelltype[matched] <- mapping_vector[meta_cells[matched]]

dir.create("niche4_start")
# (3). Calculate the cell-cell communication scores
for (cancer_type in unique(allcancer_niche4@meta.data$cancertype)) {
  print(cancer_type)
  seurat_split=eval(parse(text = paste0("subset(allcancer_niche4,subset = cancertype == \"",cancer_type,"\")")))
  Idents(seurat_split) <- "majorcelltype"
  seurat_split
  
  db <- CCC_data
  # Parse ligand-receptor pairs
  parse_lr <- function(x) {
    parts <- str_split_fixed(x, "_", 2)
    lig   <- str_trim(parts[1])
    rec   <- str_trim(parts[2])
    list(ligand = lig, receptor = rec)
  }
  
  parsed <- lapply(db$interaction, parse_lr)
  db$ligand_genes   <- sapply(parsed, `[[`, "ligand")
  db$receptor_genes <- sapply(parsed, `[[`, "receptor")
  
  # Calculate average expression
  avg_expr <- AverageExpression(seurat_split, 
                                assays = "Xenium",   
                                slot   = "data", 
                                return.seurat = FALSE)$Xenium
  avg_expr <- as.data.frame(avg_expr)
  avg_expr <- t(avg_expr) %>% as.data.frame()
  celltypes <- rownames(avg_expr)
  
  # Get expression of each L-R in each cell type
  L_expr <- matrix(0, nrow = length(celltypes), ncol = nrow(db),
                   dimnames = list(celltypes, db$interaction))
  R_expr <- matrix(0, nrow = length(celltypes), ncol = nrow(db),
                   dimnames = list(celltypes, db$interaction))
  
  for (i in 1:nrow(db)) {
    if (db$ligand_genes[[i]] %in% colnames(avg_expr) & db$receptor_genes[[i]] %in% colnames(avg_expr)) {
      L_expr[, i] <- avg_expr[,db$ligand_genes[[i]]]
      R_expr[, i] <- avg_expr[,db$receptor_genes[[i]]]
    }
  }
  
  # Calculate communication strength (source -> target)
  results <- expand.grid(
    source      = celltypes,
    target      = celltypes,
    interaction = db$interaction,
    stringsAsFactors = FALSE
  ) %>% mutate(score = 0)
  
  for (src in celltypes) {
    for (tgt in celltypes) {
      score_vec <- sqrt( L_expr[src, ] * R_expr[tgt, ] )
      idx <- results$source == src & results$target == tgt
      results$score[idx] <- score_vec
    }
  }
  
  # Save results to CSV
  output_file <- paste0("~/space/ccc/niche4_start/cellchat_", cancer_type, ".net.csv")
  
  write.csv(results, file = output_file, row.names = FALSE)
  
}

################################################################################
# 4_2. Load ccc results of STAR-T and Other-T
################################################################################

# (1). Load data
tumors = c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for (j in tumors) {
  eval(parse(text = paste0("df.net.cellchat_",j,"_niche4_major=read.csv(\"~/space/ccc/niche4_start/cellchat_",j,".net.csv\")")))
  eval(parse(text = paste0("df.net.cellchat_",j,"_niche4_major$niche=\"SN4\"")))
  eval(parse(text = paste0("df.net.cellchat_",j,"_niche4_major$cancer=\"",j,"\"")))
}

# (2). Combine data of all cancers
major_cellchat_sn4=ls()[grep("df.net.cellchat",ls())]
major_cellchat_sn4=major_cellchat_sn4[grep("niche4",major_cellchat_sn4)]
major_cellchat_sn4
df.net.cellchat_sn4=do.call(rbind,obj_list <- mget(major_cellchat_sn4, envir = .GlobalEnv))
df.net.cellchat_sn4=df.net.cellchat_sn4[order(df.net.cellchat_sn4$cancer),]
row.names(df.net.cellchat_sn4)=NULL

# (3). Set the group of STAR-T vs Other-T
df_filtered_sn4=df.net.cellchat_sn4[df.net.cellchat_sn4$niche=="SN4",]
head(df_filtered_sn4)

# (4). Standardize cell type names of STAR-T and Other-T
df_filtered_sn4$target[df_filtered_sn4$target=="cell-other"]="cell_other"
df_filtered_sn4$target[df_filtered_sn4$target=="cell-start"]="cell_start"
df_filtered_sn4$source[df_filtered_sn4$source=="cell-other"]="cell_other"
df_filtered_sn4$source[df_filtered_sn4$source=="cell-start"]="cell_start"
unique(df_filtered_sn4$source)

df_filtered_sn4$source=gsub("T-cell","T_cell",df_filtered_sn4$source)
df_filtered_sn4$source=gsub("B-cell","B_cell",df_filtered_sn4$source)
df_filtered_sn4$source=gsub("Dendritic-cell","Dendritic_cell",df_filtered_sn4$source)
df_filtered_sn4$source=gsub("Smooth-muscle","Smooth_muscle",df_filtered_sn4$source)
df_filtered_sn4$target=gsub("T-cell","T_cell",df_filtered_sn4$target)
df_filtered_sn4$target=gsub("B-cell","B_cell",df_filtered_sn4$target)
df_filtered_sn4$target=gsub("Dendritic-cell","Dendritic_cell",df_filtered_sn4$target)
df_filtered_sn4$target=gsub("Smooth-muscle","Smooth_muscle",df_filtered_sn4$target)
unique(df_filtered_sn4$source)

# (5). Append a pair column, include source-target-ligand-receptor
df_filtered_sn4$source_target=paste(df_filtered_sn4$source,df_filtered_sn4$target,sep=":")
df_filtered_sn4$pair=paste(df_filtered_sn4$source_target,df_filtered_sn4$interaction,sep=":")
colnames(df_filtered_sn4)[4]="prob"

head(df_filtered_sn4)
df_filtered_sn4$ligand=sapply(strsplit(df_filtered_sn4$interaction, split = "_"),function(x){c(x[[1]])})
df_filtered_sn4$receptor=sapply(strsplit(df_filtered_sn4$interaction, split = "_"),function(x){c(x[[2]])})

rm(list=paste0("df.net.cellchat_",tumors,"_niche4_major"))
gc()

# (7). Load single-cell reference data containing the expression proportion of each gene in each cell type
# Gene expression percentage data were derived from the public pan-cancer scRNA-seq dataset reported in "Cross-tissue multicellular coordination and its rewiring in cancer"
sc_ref_percentage=read.csv("~/space/analysis/signature/pancancer_celllineage.genes.percent.csv")
unique(sc_ref_percentage$Celltype)

# (8). define a genelist of cell prop > 5%
marker_list <- sc_ref_percentage %>%
  group_by(Celltype) %>%
  arrange(desc(Proportion)) %>%
  summarise(
    gene_gt05 = list(Gene[Proportion > 0.05]),  
    .groups = "drop"
  )

final_marker_list <- marker_list %>%
  { setNames(.$gene_gt05, paste0(.$Celltype, "_gt05%")) }

allgenes_5k=read.csv("~/space/spatial_5k_genes.csv")

selected_marker_genes <- list(
  Endothelial_gt05 = unlist(final_marker_list[["Endothelial_gt05%"]]),
  Fibroblast_gt05 = unlist(final_marker_list[["Fibroblast_gt05%"]]),
  Macrophage_gt05 = unlist(final_marker_list[["Macrophage_gt05%"]]),     
  T_cell_gt05 = unlist(final_marker_list[["T_cell_gt05%"]]),
  Pericyte_gt05 = unlist(final_marker_list[["Pericyte_gt05%"]]),
  Plasma_gt05 = unlist(final_marker_list[["Plasma_gt05%"]]),
  B_cell_gt05 = unlist(final_marker_list[["B_cell_gt05%"]]),
  Smooth_muscle_gt05 = unlist(final_marker_list[["Smooth_muscle_gt05%"]]),
  Mast_gt05 = unlist(final_marker_list[["Mast_gt05%"]]),
  Dendritic_cell_gt05 = unlist(final_marker_list[["Dendritic_cell_gt05%"]]),
  
  Malignant_gt05 = allgenes_5k$gene, 
  cell_other_gt05 = unlist(final_marker_list[["T_cell_gt05%"]]),
  cell_start_gt05 = unlist(final_marker_list[["T_cell_gt05%"]])
  
)

# (9). Filtering interactions with at least 5% expression in corresponding cell types of single cell reference dataset
setDT(df_filtered_sn4)
marker_genes <- lapply(selected_marker_genes, as.character)
dt <- as.data.table(df_filtered_sn4)
dt=dt[dt$prob!=0,]
dt=dt[!(dt$source %in% c("Oligodendrocyte","Keratinocyte","Neural","Neutrophil","Monocyte")) & !(dt$target %in% c("Oligodendrocyte","Keratinocyte","Neural","Neutrophil","Monocyte")),]
unique(dt$source)

df_filtered_sn4_prop0.05 <- dt[
  # Ligand gene is expressed >5% in the source cell
  , ligand_expressed := mapply(function(genes, cell) any(genes %in% marker_genes[[paste0(cell, "_gt05")]]),
                               ligand, source)
][
  # Receptor gene is expressed >5% in the target cell
  , receptor_expressed := mapply(function(genes, cell) any(genes %in% marker_genes[[paste0(cell, "_gt05")]]),
                                 receptor, target)
][
  # Final filtering
  ligand_expressed == TRUE & receptor_expressed == TRUE
][
  # Clean up temporary columns
  , `:=`(ligand = NULL, receptor = NULL,
         ligand_expressed = NULL, receptor_expressed = NULL)
]

cat("Remaining after filtering:", nrow(df_filtered_sn4_prop0.05), "high-confidence interactions (both L and R are highly expressed in corresponding cells)\n")

# (10). Export
write.csv(df_filtered_sn4_prop0.05,file="~/space/ccc/4_2_allccc_pairs_start_refprop0.05.csv")


################################################################################
# 4_3. Paired test to identify up and down-regulated ligand-receptor pairs in STAR-T 
################################################################################

# (1). Transfer long data into wide matrix for paired test
# Retain only t cells act as source or target
df_filtered_sn4_prop0.05_t = df_filtered_sn4_prop0.05[df_filtered_sn4_prop0.05$source %in% c("cell_other","cell_start") | df_filtered_sn4_prop0.05$target %in% c("cell_other","cell_start"),]
df_filtered_sn4_prop0.05_t=df_filtered_sn4_prop0.05_t[-intersect(grep("cell_", df_filtered_sn4_prop0.05_t$source),grep("cell_", df_filtered_sn4_prop0.05_t$target)) ,]
# Use Niche to record whether it is cell_start or cell_other
df_filtered_sn4_prop0.05_t[, niche := ifelse(source %in% c("cell_start", "cell_other"), source, target)]

# Set the group of T cells
df_filtered_sn4_prop0.05_t[, c("source", "target") := lapply(.SD, function(x) {
  x[x %in% c("cell_start", "cell_other")] <- "T_cell"
  x
}), .SDcols = c("source", "target")]
head(df_filtered_sn4_prop0.05_t)
# Replace STAR-T/Other-T in pair with T cell
df_filtered_sn4_prop0.05_t[, pair := gsub("cell_start|cell_other", "T_cell", pair)]
head(df_filtered_sn4_prop0.05_t)

# Covert into wide format
df <- df_filtered_sn4_prop0.05_t
df$pair <- as.factor(df$pair)
dt <- as.data.table(df)
dt$niche_cancer=paste(dt$niche,dt$cancer,sep="_")
system.time({
  df_wide_dt <- dcast(
    dt,
    niche + cancer + niche_cancer ~ pair,
    value.var = "prob",
    fill = 0, 
    fun.aggregate = sum
  )
})

df_wide_dt=as.data.frame(df_wide_dt)
ccc_sn4_tcells_wide <-  df_wide_dt %>% arrange(niche, cancer)

# (2). Export the ccc matrix
write.csv(ccc_sn4_tcells_wide,file="~/space/ccc/4_3_ccc_start_wide.csv")

# (3). Perform paired test
df <- ccc_sn4_tcells_wide
all_SN <- sort(unique(df$niche))

df_long <- df %>%
  pivot_longer(-c(niche, cancer, niche_cancer), names_to = "pair", values_to = "prob") %>%
  as.data.table()

setDT(df_long)  # Ensure it is a data.table
head(df_long)

result <- rbindlist(lapply(all_SN, function(focal) {
  message("Calculating: ", focal, " vs other niches")
  
  # First calculate focal and other values by pair + cancer
  tmp <- df_long[, .(
    prob_focal = ifelse(any(niche == focal), prob[niche == focal], NA),
    prob_other = mean(prob[niche != focal])   # Average of other niches in the same cancer
  ), by = .(pair, cancer)]
  
  # Remove NA values (some cancers may not have this focal niche)
  tmp <- tmp[!is.na(prob_focal)]
  
  # Summarize by pair and perform paired test
  tmp[, .(
    focal_SN     = focal,
    n_cancer     = .N,                                 # Number of comparable cancers
    mean_focal   = mean(prob_focal),
    mean_other   = mean(prob_other),
    fold_change  = mean(prob_focal) / (mean(prob_other) + 1e-8),
    mean_diff    = mean(prob_focal) - mean(prob_other),
    
    p_wilcox     = if(.N >= 3) wilcox.test(prob_focal, prob_other, paired = TRUE)$p.value else NA_real_
  ), by = pair]
}))

# (4). Append significance stars and sort the results
result <- result %>%
  mutate(sig = case_when(
    p_wilcox < 0.001 ~ "***",
    p_wilcox < 0.01  ~ "**",
    p_wilcox < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  arrange(focal_SN, p_wilcox)

sum(result$p_wilcox<0.05)
result=result[!is.na(result$p_wilcox),]
ccc_sn4_tcells_pairedtest=result

# (5). Extract significant results
ccc_sn4_tcells_pairedtest_sig=ccc_sn4_tcells_pairedtest[ccc_sn4_tcells_pairedtest$p_wilcox<0.05,]
ccc_sn4_tcells_pairedtest_sig=ccc_sn4_tcells_pairedtest_sig[ccc_sn4_tcells_pairedtest_sig$focal_SN=="cell_start",]

# (6). Save the results
write.csv(ccc_sn4_tcells_pairedtest_sig,file="~/space/ccc/4_3_ccc_start_pairedtest_sig.csv")


