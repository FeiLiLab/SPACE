library(Seurat)
library(dplyr)
library(BiocParallel)
library(SeuratWrappers)
library(future)
options(future.globals.maxSize = 300 * 1024^3)


################# 1. Identify spatially enriched genes (SEGs) for each cell type across cancer types and spatial niches
#################
#################
#################
#################
#################


################################################################################
# 1_1. Load data
################################################################################

seurat_niches<-readRDS("~/space/seurat_niche/metadata_merge.rds")  # Merged meta data of all cancers was generated in Step02, line 139

cancer_list <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for(cancer in cancer_list){
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS files of each cancer were generated in Step02, line 135
  assign(paste0("obj_", cancer), obj)
  message("Processed: ", cancer)
}

################################################################################
# 1_2. Identify spatially enriched genes (SEGs)
################################################################################
plan("multisession",workers=16)
deg_all_target<-data.frame()
deg_all_other<-data.frame()
deg_all_overlap<-data.frame()
majorcelltype_list<-setdiff(names(table(seurat_niches$majorcelltype)),c("Keratinocyte","Oligodendrocyte"))
for (majorcelltype_select in majorcelltype_list){
  seurat_niches_target<-seurat_niches[which(seurat_niches$majorcelltype==majorcelltype_select),]
  # Filter cancer types with low cell number
  cancer_types <- seurat_niches_target %>%
    group_by(cancertype, niche) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(cancertype) %>%
    summarise(all_above_3 = all(count > 2)) %>%
    filter(all_above_3) %>%
    pull(cancertype)
  # Identify DEGs
  for (cancer in cancer_types) {
    obj <- get(paste0("obj_", cancer))
    merged_metadata <- obj@meta.data %>%
      tibble::rownames_to_column("cell_id") %>%
      left_join(tibble::rownames_to_column(seurat_niches_target["niche"], "cell_id"),by = "cell_id") %>%
      tibble::column_to_rownames("cell_id")
    obj$niche_target <- ifelse(is.na(merged_metadata$niche.y),"Other",merged_metadata$niche.y)
    
    Idents(obj)<-obj$niche_target
    niche_list<-setdiff(names(table(obj$niche_target)),"Other")
    for (niche_select in niche_list){
      # Identify the DEGs between niche_select and other niches within the same cell type
      deg_all_target_1<-FindMarkers(obj, ident.1 = niche_select,ident.2 = setdiff(names(table(obj$niche_target)),c(niche_select,"Other")))
      deg_all_target_1$cancertype<-rep(cancer,nrow(deg_all_target_1))
      deg_all_target_1$niche<-rep(niche_select,nrow(deg_all_target_1))
      deg_all_target_1$gene<-row.names(deg_all_target_1)
      deg_all_target_1$majorcelltype<-rep(majorcelltype_select,nrow(deg_all_target_1))
      # Identify the DEGs between niche_select and other different cell types
      deg_all_other_1<-FindMarkers(obj, ident.1 = niche_select,ident.2 = c("Other"))
      deg_all_other_1$cancertype<-rep(cancer,nrow(deg_all_other_1))
      deg_all_other_1$niche<-rep(niche_select,nrow(deg_all_other_1))
      deg_all_other_1$gene<-row.names(deg_all_other_1)
      deg_all_other_1$majorcelltype<-rep(majorcelltype_select,nrow(deg_all_other_1))
      
      # Identify the DEGs overlapped between the above two comparisons (focusing on up-regulated DEGs only)
      deg_all_overlap_1<-deg_all_target_1[intersect(row.names(deg_all_target_1[which(deg_all_target_1$avg_log2FC>0&deg_all_target_1$p_val_adj<0.05),]),
                                                    row.names(deg_all_other_1[which(deg_all_other_1$avg_log2FC>0&deg_all_other_1$p_val_adj<0.05),])),]
      # Combine results from different cancer types
      deg_all_target<-rbind(deg_all_target,deg_all_target_1)
      deg_all_other<-rbind(deg_all_other,deg_all_other_1)
      deg_all_overlap<-rbind(deg_all_overlap,deg_all_overlap_1)
      message(paste("Processed", cancer,"-",majorcelltype_select,"-",niche_select,":", nrow(deg_all_overlap_1), "DEGs."))
    }
  }
}
dir.create("~/space/analysis/domain/deg/",recursive = T)
write.csv(deg_all_target,"~/space/analysis/domain/deg/1_2_deg_all_target.csv")
write.csv(deg_all_other,"~/space/analysis/domain/deg/1_2_deg_all_other.csv")
write.csv(deg_all_overlap,"~/space/analysis/domain/deg/1_2_deg_all_overlap.csv")


################################################################################
# 1_3. Retain candidate genes that were consistently detected (Percentage > 5%) in the corresponding cell type in the public pan-cancer scRNA-seq dataset
################################################################################

# Load data 
# Gene expression percentage data were derived from the public pan-cancer scRNA-seq dataset reported in "Cross-tissue multicellular coordination and its rewiring in cancer"
pancancer_gene<-read.csv("~/space/analysis/signature/pancancer_celllineage.genes.percent.csv")
unique(pancancer_gene$Celltype)

# Identify overlapped cell types
celltype_overlap<-intersect(names(table(deg_all_overlap$majorcelltype)),names(table(pancancer_gene$Celltype)))
celltype_analysis<-setdiff(celltype_overlap,c("Monocyte","Neutrophil"))

# Perform gene filtering
deg_final<-data.frame()
for (celltype in celltype_analysis){
  pancancer_gene_select<-pancancer_gene[which(pancancer_gene$Celltype==celltype&pancancer_gene$Proportion>0.05),]
  gene_expressed<-pancancer_gene_select$Gene 
  deg_final_1<-deg_all_overlap[which(deg_all_overlap$majorcelltype==celltype&deg_all_overlap$gene %in% gene_expressed),]
  deg_final<-rbind(deg_final,deg_final_1)
}
write.csv(deg_final,"~/space/analysis/domain/deg/1_3_deg_final.csv")


################################################################################
# 1_4. Identify cancer-shared SEGs
################################################################################

deg_final<-read.csv("~/space/analysis/domain/deg/1_3_deg_final.csv")
deg_final$cluster<-deg_final$niche
deg_final$celltype<-deg_final$majorcelltype

gene_cancertype_count <- deg_final %>%
  group_by(cluster, celltype, gene) %>%
  summarise(cancertype_count = n_distinct(cancertype), .groups = 'drop')

gene_at_least_3_cancertypes <- gene_cancertype_count %>%
  filter(cancertype_count >= 3)

deg_merge_crosscancer <- deg_final %>%
  inner_join(gene_at_least_3_cancertypes, 
             by = c("cluster", "celltype", "gene"))
write.csv(deg_merge_crosscancer,file="~/space/analysis/domain/deg/1_4_deg_merge_crosscancer.csv")




################# 2. Pathway enrichment of SEGs
#################
#################
#################
#################
#################

################################################################################
# 2_1. Load SEGs
################################################################################

deg_merge_crosscancer<-read.csv("~/space/analysis/domain/deg/1_4_deg_merge_crosscancer.csv",row.names = 1)

################################################################################
# 2_2. Perform reactome pathway enrichment
################################################################################
library(clusterProfiler)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(biomaRt)
library(ReactomePA)

data_reactome_merge<-data.frame()
celltype_list<-names(table(deg_merge_crosscancer$majorcelltype))
for (celltype in celltype_list){

  deg_select<-deg_merge_crosscancer[which(deg_merge_crosscancer$majorcelltype==celltype),]
  niche_list<-names(table(deg_select$niche))
  for (niche in niche_list){
    tryCatch({
      gene_select<-unique(deg_select[which(deg_select$niche==niche),]$gene)
      geneid_select <- bitr(gene_select, fromType = "SYMBOL",toType = c("ENSEMBL", "SYMBOL","ENTREZID"),OrgDb = org.Hs.eg.db) 
      reactome <- enrichPathway(gene=geneid_select$ENTREZID,pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable=TRUE)
      data_reactome<-data.frame(reactome)
      data_reactome$neglogpvalue<-(-log(data_reactome$pvalue,10))
      data_reactome$neglogqvalue<-(-log(data_reactome$qvalue,10))
      data_reactome_unique <-as.data.frame(data_reactome %>% group_by(Description) %>% filter(row_number() == which.min(qvalue)) %>% ungroup())
      data_reactome_unique$Description<-factor(data_reactome_unique$Description,levels=rev(data_reactome_unique$Description))
      data_reactome_unique$celltype<-rep(celltype,nrow(data_reactome_unique))
      data_reactome_unique$niche<-rep(niche,nrow(data_reactome_unique))
      data_reactome_merge<-rbind(data_reactome_unique,data_reactome_merge)
    }, error = function(e) {
      message("Process niche ", niche, " Error: ", e$message)
    })
    message(paste("Processed", celltype,"-",niche,":", nrow(data_reactome_unique), "Pathways."))
  }
}
write.csv(data_reactome_merge,"~/space/analysis/domain/deg/2_2_data_reactome_merge.csv")
