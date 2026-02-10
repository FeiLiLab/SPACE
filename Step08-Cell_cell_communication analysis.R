library(Seurat)
library(Matrix)
library(tidyr)
library(stringr)
library(dplyr)
library(tibble)
library(data.table)
library(tidyverse)
options(future.globals.maxSize = 300 * 1024^3)



#################1. Perform niche-specific analysis of cell-cell communications
#################
#################
#################
#################
#################
#################

################################################################################
# 1_1. Prepare data for cell-cell communications analysis
################################################################################

dir.create("~/space/ccc/")
setwd("~/space/ccc/")

cancer_list <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for(cancer in cancer_list){
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS files of each cancer were generated in Step02, line 135
  assign(paste0("obj_", cancer), obj)
  message("Processed: ", cancer)
}

# Merge data of all cancers into a combined Seurat object
cancerobj_list=mget(paste0("obj_", cancer_list))
merged_obj=merge(x=cancerobj_list[[1]],y=cancerobj_list[-1])

# Remove temporary files to reduce memory usage 
rm(cancerobj_list)
gc()
for(cancer in cancer_list){
  eval(parse(text = paste0("rm(obj_",cancer,")")))
}
gc()

# Set the image part of Seurat object to empty to reduce memory usage 
merged_obj@images=list()

# Extract data of each spatial niche and export
niche_numbers <- 1:5
for (i in niche_numbers) {
  current_seurat_name <- paste0("allcancer_niche", i)
  
  allcancer_currentniche <- subset(merged_obj, subset = niche == i)
  
  DefaultAssay(allcancer_currentniche)
  allcancer_currentniche <- JoinLayers(allcancer_currentniche)
  
  assign(paste0(current_seurat_name,""), allcancer_niche)
  
  saveRDS(allcancer_currentniche,file=paste0(current_seurat_name,".rds"))
  cat(current_seurat_name, "saved successfully!\n\n")
}


################################################################################
# 1_2. Calculate the cell-cell communication scores
################################################################################

niche_numbers <- 1:5
# Load CCC ligand-receptor data
# Human ligand-receptor interactions from CellChat (version 2.2.0), CellPhoneDB (version 5.0.0), and iTALK were merged, and only one-to-one L-R pairs were retained.
CCC_data <- read.csv("~/space/ccc/cellchatdb_cpdb_italk.csv")

# Loop through each niche
for (n in niche_numbers) {
  # Load the corresponding Seurat object
  niche_name <- paste0("allcancer_niche", n)
  seurat_obj_niche <- get(niche_name)
  # Create output directory
  output_dir <- paste0("./niche", n, "/")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each cancer type for current niche
  for (cancer_type in unique(seurat_obj_niche@meta.data$cancertype)) {
    # Subset Seurat object by cancer type
    print(cancer_type)
    seurat_split <- subset(seurat_obj_niche, subset = cancertype == cancer_type)
    print(dims(seurat_split))
    Idents(seurat_split) <- "majorcelltype"
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
    output_file <- paste0(output_dir, "cellchat_", cancer_type, ".net.csv")
    write.csv(results, file = output_file, row.names = FALSE)
  }
}



#################2. Pre-process results of cell-cell communication analysis
#################
#################
#################
#################
#################
#################

################################################################################
# 2_1. Load and pre-process of cell-cell communication results
################################################################################

# (1). Set parameters
setwd("~/space/ccc/")
tumors = c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")

# (2). Load data of SN1 to SN4
for (i in 1:4) {
  for (j in tumors) {
    eval(parse(text = paste0("df.net.cellchat_",j,"_niche",i,"_major=read.csv(\"~/space/ccc/niche",i,"/cellchat_",j,".net.csv\")")))
    eval(parse(text = paste0("df.net.cellchat_",j,"_niche",i,"_major$niche=\"SN",i,"\"")))
    eval(parse(text = paste0("df.net.cellchat_",j,"_niche",i,"_major$cancer=\"",j,"\"")))
  }
}

# (3). Load data of SN5, as some cancer don't have SN5
dir("~/space/ccc/niche5/")
tumors_niche5=dir("~/space/ccc/niche5/")
tumors_niche5=tumors_niche5[grep("csv",tumors_niche5)]
tumors_niche5=gsub("cellchat_","",tumors_niche5)
tumors_niche5=gsub(".net.csv","",tumors_niche5)
tumors_niche5

for (i in 5) {
  for (j in tumors_niche5) {
    eval(parse(text = paste0("df.net.cellchat_",j,"_niche5_major=read.csv(\"~/space/ccc/niche5/cellchat_",j,".net.csv\")")))
    eval(parse(text = paste0("df.net.cellchat_",j,"_niche5_major$niche=\"SN5\"")))
    eval(parse(text = paste0("df.net.cellchat_",j,"_niche5_major$cancer=\"",j,"\"")))
  }
}

# (4). Combine all cancer and niches
major_cellchat_all=ls()[grep("df.net.cellchat",ls())]
major_cellchat_all=major_cellchat_all[grep("niche",major_cellchat_all)]
major_cellchat_all
df.net.cellchat_all_major=do.call(rbind,obj_list <- mget(major_cellchat_all, envir = .GlobalEnv))

df.net.cellchat_all_major=df.net.cellchat_all_major[order(df.net.cellchat_all_major$cancer),]
row.names(df.net.cellchat_all_major)=NULL
head(df.net.cellchat_all_major)

# (5). Append a column of pair, which including source-target-ligand-receptor
# Set the separator of "-" to "_", as Seurat will auto-change "_" to  "-"
unique(df.net.cellchat_all_major$source)
df.net.cellchat_all_major$source=gsub("T-cell","T_cell",df.net.cellchat_all_major$source)
df.net.cellchat_all_major$source=gsub("B-cell","B_cell",df.net.cellchat_all_major$source)
df.net.cellchat_all_major$source=gsub("Dendritic-cell","Dendritic_cell",df.net.cellchat_all_major$source)
df.net.cellchat_all_major$source=gsub("Smooth-muscle","Smooth_muscle",df.net.cellchat_all_major$source)
df.net.cellchat_all_major$target=gsub("T-cell","T_cell",df.net.cellchat_all_major$target)
df.net.cellchat_all_major$target=gsub("B-cell","B_cell",df.net.cellchat_all_major$target)
df.net.cellchat_all_major$target=gsub("Dendritic-cell","Dendritic_cell",df.net.cellchat_all_major$target)
df.net.cellchat_all_major$target=gsub("Smooth-muscle","Smooth_muscle",df.net.cellchat_all_major$target)

df.net.cellchat_all_major$pair=paste(df.net.cellchat_all_major$source,df.net.cellchat_all_major$target,sep=":")
df.net.cellchat_all_major$pair=paste(df.net.cellchat_all_major$pair,df.net.cellchat_all_major$interaction,sep=":")

# (6). Extract data of 10 cell types
sum(c("T_cell","Plasma","B_cell","Macrophage","Fibroblast","Smooth_muscle","Dendritic_cell","Pericyte","Endothelial","Mast") %in% unique(df.net.cellchat_all_major$source))
df.net.cellchat_all_major=df.net.cellchat_all_major[df.net.cellchat_all_major$source %in% c("T_cell","Plasma","B_cell","Macrophage","Fibroblast","Smooth_muscle","Dendritic_cell","Pericyte","Endothelial","Mast"),]
df.net.cellchat_all_major=df.net.cellchat_all_major[df.net.cellchat_all_major$target %in% c("T_cell","Plasma","B_cell","Macrophage","Fibroblast","Smooth_muscle","Dendritic_cell","Pericyte","Endothelial","Mast"),]


# (7). Standardize the ligand and receptor names of data
df.net.cellchat_all_major$ligand=sapply(strsplit(df.net.cellchat_all_major$interaction, split = "_"),function(x){c(x[[1]])})
df.net.cellchat_all_major$receptor=sapply(strsplit(df.net.cellchat_all_major$interaction, split = "_"),function(x){c(x[[2]])})
colnames(df.net.cellchat_all_major)[4]="prob"

# (8). Remove temporary files
for (i in 1:5) {
  for (j in tumors) {
    eval(parse(text = paste0("rm(df.net.cellchat_",j,"_niche",i,"_major)")))
  }
}
rm(obj_list)
head(df.net.cellchat_all_major)


################################################################################
# 2_2. Data clean of ccc (cell-cell communication) results
################################################################################

# (1). Load single-cell reference data containing the expression proportion of each gene in each cell type
# Gene expression percentage data were derived from the public pan-cancer scRNA-seq dataset reported in "Cross-tissue multicellular coordination and its rewiring in cancer"
sc_ref_percentage=read.csv("~/space/analysis/signature/pancancer_celllineage.genes.percent.csv")
unique(sc_ref_percentage$Celltype)

# (3). Define a genelist of cell prop > 5%
marker_list <- sc_ref_percentage %>%
  group_by(Celltype) %>%
  arrange(desc(Proportion)) %>%
  summarise(
    gene_gt05 = list(Gene[Proportion > 0.05]),  
    .groups = "drop"
  )

final_marker_list <- marker_list %>%
  { setNames(.$gene_gt05, paste0(.$Celltype, "_gt05%")) }


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
  Dendritic_cell_gt05 = unlist(final_marker_list[["Dendritic_cell_gt05%"]])
)

# (4). Filtering interactions with at least 5% expression in corresponding cell types of single cell reference dataset
setDT(df.net.cellchat_all_major)

marker_genes <- lapply(selected_marker_genes, as.character)

dt <- as.data.table(df.net.cellchat_all_major)
dt=dt[dt$prob!=0,]

head(dt)
df_filtered_ref_prop0.05 <- dt[
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

cat("Remaining after filtering:", nrow(df_filtered_ref_prop0.05), "high-confidence interactions (both L and R are highly expressed in corresponding cells)\n")

# (5). Export
write.csv(df_filtered_ref_prop0.05,file="~/space/ccc/2_2_allccc_pairs_refprop0.05.csv")


################################################################################
# 2_3. Paired test to identify SN specify up-regulated ccc pairs
################################################################################

# (1). Transfer long data into wide matrix for paired test
df <- df_filtered_ref_prop0.05
df$pair <- as.factor(df$pair)

dt <- as.data.table(df)
dt$niche_cancer=paste(dt$niche,dt$cancer,sep="_")

system.time({
  df_wide_dt <- dcast(
    dt,
    niche + cancer + niche_cancer ~ pair,   # Columns on the left are retained, column on the right is spread wide
    value.var = "prob",
    fill = 0,                               # Fill missing values with 0
    fun.aggregate = sum                     # Sum values for duplicates (duplicates are generally not present)
  )
})

df_wide_dt=as.data.frame(df_wide_dt)
df_ref_prop0.05_allpairs_wide <-  df_wide_dt %>% arrange(niche, cancer)

# (2). Fill the ccc of missed cancers in SN5
setdiff(df_ref_prop0.05_allpairs_wide$cancer[df_ref_prop0.05_allpairs_wide$niche=="SN1"],df_ref_prop0.05_allpairs_wide$cancer[df_ref_prop0.05_allpairs_wide$niche=="SN5"])
#[1] "brca" "gbm"  "lihc"
df_ref_prop0.05_allpairs_wide=rbind(df_ref_prop0.05_allpairs_wide,df_ref_prop0.05_allpairs_wide[1:3,])
tail(df_ref_prop0.05_allpairs_wide[,1:5])
df_ref_prop0.05_allpairs_wide[(nrow(df_ref_prop0.05_allpairs_wide)-2):(nrow(df_ref_prop0.05_allpairs_wide)),1]="SN5"
df_ref_prop0.05_allpairs_wide[(nrow(df_ref_prop0.05_allpairs_wide)-2):(nrow(df_ref_prop0.05_allpairs_wide)),2]=c("brca","gbm","lihc")
df_ref_prop0.05_allpairs_wide[(nrow(df_ref_prop0.05_allpairs_wide)-2):(nrow(df_ref_prop0.05_allpairs_wide)),3]=paste(df_ref_prop0.05_allpairs_wide[(nrow(df_ref_prop0.05_allpairs_wide)-2):(nrow(df_ref_prop0.05_allpairs_wide)),1],df_ref_prop0.05_allpairs_wide[(nrow(df_ref_prop0.05_allpairs_wide)-2):(nrow(df_ref_prop0.05_allpairs_wide)),2],sep="_")
tail(df_ref_prop0.05_allpairs_wide[,1:5])

# (3). Set NA to 0
df_ref_prop0.05_allpairs_wide[(nrow(df_ref_prop0.05_allpairs_wide)-2):(nrow(df_ref_prop0.05_allpairs_wide)),4:ncol(df_ref_prop0.05_allpairs_wide)]=0
# (4). Sort
df_ref_prop0.05_allpairs_wide <-  df_ref_prop0.05_allpairs_wide %>% arrange(niche, cancer)
tail(df_ref_prop0.05_allpairs_wide[,1:5])

# (5). Export
write.csv(df_ref_prop0.05_allpairs_wide,file="~/space/ccc/2_3_allccc_pairs_refprop0.05_wide.csv")

# (6). Paired test
df <- df_ref_prop0.05_allpairs_wide
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
    prob_focal = ifelse(any(niche == focal), prob[niche == focal], as.numeric(NA)),
    prob_other = ifelse(any(niche != focal), mean(prob[niche != focal]), as.numeric(NA))   # Average of other niches in the same cancer
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
    
    # Now it's safe to perform the paired test!
    p_wilcox     = if(.N >= 3) wilcox.test(prob_focal, prob_other, paired = TRUE)$p.value else NA_real_
  ), by = pair]
}))

# (7). Append significance stars and sort the results
result <- result %>%
  mutate(sig = case_when(
    p_wilcox < 0.001 ~ "***",
    p_wilcox < 0.01  ~ "**",
    p_wilcox < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  arrange(focal_SN, p_wilcox)

result=result[!is.na(result$p_wilcox),]

# (8). Extract significant up-regulated pairs in each SN
ccc_prop0.05_pairedtest = result
ccc_prop0.05_pairedtest_sig=ccc_prop0.05_pairedtest[ccc_prop0.05_pairedtest$p_wilcox<0.05 & ccc_prop0.05_pairedtest$mean_diff>0,]

ccc_prop0.05_pairedtest_sig$source=sapply(strsplit(ccc_prop0.05_pairedtest_sig$pair, split = ":"),function(x){c(x[[1]])})
ccc_prop0.05_pairedtest_sig$target=sapply(strsplit(ccc_prop0.05_pairedtest_sig$pair, split = ":"),function(x){c(x[[2]])})
ccc_prop0.05_pairedtest_sig$LR=sapply(strsplit(ccc_prop0.05_pairedtest_sig$pair, split = ":"),function(x){c(x[[3]])})
ccc_prop0.05_pairedtest_sig$source_target = paste(ccc_prop0.05_pairedtest_sig$source,ccc_prop0.05_pairedtest_sig$target,sep=":")

# (9). Export the final result
write.csv(ccc_prop0.05_pairedtest_sig,file="~/space/ccc/2_3_ccc_prop0.05_pairedtest_sig.csv")


