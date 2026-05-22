options(future.globals.maxSize = 300 * 1024^3)
library(Seurat)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(BiocParallel)
library(harmony)
library(patchwork)
library(tidyverse)
library(paletteer)
library(hdf5r)
library(presto)
library(future)
library(ggpubr)
library(DESeq2)
library(pheatmap)
library(ggrepel)
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


################# 1. Analysis of T cells in Gain-of-function-based Perturb-seq of 36 STAR-T signature genes
################# 
#################
#################
#################
#################

################################################################################
# 1. Data Load, preprocess and QC
################################################################################

library(Seurat)
library(Matrix)

# (1). Load the Seurat object of BD Rhapsody output, B16OT1_5
for (i in 5:8) {
  assign(paste0("B16OT1_",i), readRDS(paste0("~space/genes_oe/output/B16OT1_",i,"/Taglib_B16OT1_",i,"_Seurat.rds")))
}
table(B16OT1_1$Sample_Tag)


# (2). Multiplet and Undetermined filter
samples <- paste0("B16OT1_",5:8)
for (s in samples) {
  obj <- get(s)
  cells_to_keep <- colnames(obj)[!(obj@meta.data$Sample_Tag %in% c("Multiplet", "Undetermined"))]
  obj_sub <- subset(obj, cells = cells_to_keep)
  assign(s, obj_sub)
}

# (3). Filter genes in less than 3 cells
samples <- paste0("B16OT1_",5:8)
for (s in samples) {
  
  old_obj <- get(s)
  counts <- GetAssayData(old_obj, assay = "RNA", slot = "counts")
  metadata <- old_obj@meta.data
  new_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = s,min.cells = 3
  )
  assign(s, new_obj)
}

# (4). Set the cell tag to over-expressed gene names
samples <- paste0("B16OT1_",5:8)
for (s in samples) {
  obj <- get(s)
  obj$Sample_Name=as.character(obj$Sample_Name)
  assign(s, obj)
}

rename_map <- c(
  "SampleTag01_mm" = "Slamf7",
  "SampleTag02_mm" = "Adgre5",
  "SampleTag03_mm" = "Apobec3",
  "SampleTag04_mm" = "Itgae",
  "SampleTag05_mm" = "Themis",
  "SampleTag06_mm" = "Zfp683",
  "SampleTag07_mm" = "Ccr5",
  "SampleTag08_mm" = "Eomes",
  "SampleTag09_mm" = "Pdia3",
  "SampleTag10_mm" = "EV",
  "SampleTag11_mm" = "Luc"
)
B16OT1_5$Sample_Name <- recode(B16OT1_5$Sample_Name, !!!rename_map)
B16OT1_5$Sample_Name <- factor(B16OT1_5$Sample_Name, levels = as.vector(rename_map))
unique(B16OT1_5$Sample_Name)

rename_map <- c(
  "SampleTag01_mm" = "App",
  "SampleTag02_mm" = "Cxcr6",
  "SampleTag03_mm" = "Grn",
  "SampleTag04_mm" = "Itgb7",
  "SampleTag05_mm" = "S1pr5",
  "SampleTag06_mm" = "Abi3",
  "SampleTag07_mm" = "Actn4",
  "SampleTag08_mm" = "Adrb2",
  "SampleTag09_mm" = "Apmap",
  "SampleTag10_mm" = "EV",
  "SampleTag11_mm" = "Luc"
)
B16OT1_6$Sample_Name <- recode(B16OT1_6$Sample_Name, !!!rename_map)
B16OT1_6$Sample_Name <- factor(B16OT1_6$Sample_Name, levels = as.vector(rename_map))
unique(B16OT1_6$Sample_Name)

rename_map <- c(
  "SampleTag01_mm" = "Cd300a",
  "SampleTag02_mm" = "Gpr174",
  "SampleTag03_mm" = "Ifitm1",
  "SampleTag04_mm" = "Lipa",
  "SampleTag05_mm" = "Nucb2",
  "SampleTag06_mm" = "Ptgdr",
  "SampleTag07_mm" = "Ptpn22",
  "SampleTag08_mm" = "Cd38",
  "SampleTag09_mm" = "Cx3cr1",
  "SampleTag10_mm" = "EV",
  "SampleTag11_mm" = "Luc"
)
B16OT1_7$Sample_Name <- recode(B16OT1_7$Sample_Name, !!!rename_map)
B16OT1_7$Sample_Name <- factor(B16OT1_7$Sample_Name, levels = as.vector(rename_map))
unique(B16OT1_7$Sample_Name)

rename_map <- c(
  "SampleTag01_mm" = "Lpcat1",
  "SampleTag02_mm" = "Mapk1",
  "SampleTag03_mm" = "Nfatc3",
  "SampleTag04_mm" = "P4hb",
  "SampleTag05_mm" = "Ptger4",
  "SampleTag06_mm" = "Rnf167",
  "SampleTag07_mm" = "Nherf1",
  "SampleTag08_mm" = "Tbx21",
  "SampleTag12_mm" = "Tgfbr1",
  "SampleTag10_mm" = "EV",
  "SampleTag11_mm" = "Luc"
)
B16OT1_8$Sample_Name <- recode(B16OT1_8$Sample_Name, !!!rename_map)
B16OT1_8$Sample_Name <- factor(B16OT1_8$Sample_Name, levels = as.vector(rename_map))
unique(B16OT1_8$Sample_Name)

# (5). Merge into a combined seurat object
B16OT1 <- merge(
  x = B16OT1_5, 
  y = c(B16OT1_6, B16OT1_7, B16OT1_8), 
  add.cell.ids = c("S5", "S6", "S7", "S8"), 
  project = "B16OT1"
)

# (6). Add mt gene percent
B16OT1[["percent.mt"]] <- PercentageFeatureSet(B16OT1, pattern = "^mt-")

# (7). Perform QC 
pdf(file="~/space/genes_oe/1_1_B16OT1_vlnplot_beforeqc.pdf",width=12,height=4)
VlnPlot(B16OT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "Sample_Name",cols = colors_sample)
dev.off()

B16OT1<-subset(B16OT1,subset = nFeature_RNA > 200 & nCount_RNA < 40000 & percent.mt < 10)

pdf(file="~/space/genes_oe/1_1_B16OT1_vlnplot_afterqc.pdf",width=12,height=4)
VlnPlot(B16OT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "Sample_Name",cols = colors_sample)
dev.off()


################################################################################
# 2. Dimensionality reduction and clustering
################################################################################

B16OT1 <- NormalizeData(B16OT1)
B16OT1 <- FindVariableFeatures(B16OT1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(B16OT1)
B16OT1 <- ScaleData(B16OT1, features = all.genes)
B16OT1 <- RunPCA(B16OT1, features = VariableFeatures(object = B16OT1))
ElbowPlot(B16OT1)
B16OT1 <- FindNeighbors(B16OT1, dims = 1:20)
B16OT1 <- FindClusters(B16OT1, resolution = 0.1)
B16OT1 <- RunUMAP(B16OT1, dims = 1:20)

B16OT1 <- JoinLayers(B16OT1)
B16OT1.markers <- FindAllMarkers(B16OT1, only.pos = TRUE)
write.csv(B16OT1.markers,file = '~/space/genes_oe/2_1.B16OT1.markers.csv')

# Dimplot
pdf(file="~/space/genes_oe/2_1_dimplot.B16OT1.pdf",width = 12,height = 12)
DimPlot(B16OT1, reduction = "umap",label = T,cols = mypal_1)
dev.off()

pdf(file="~/space/genes_oe/2_1_dimplot.B16OT1_samples.pdf",width = 12,height = 12)
DimPlot(B16OT1, reduction = "umap",label = T,cols = colors_sample,group.by = "Sample_Name")
dev.off()

# Filter of Non-T cells
B16OT1 <- subset(B16OT1,idents = c(0,1))


################################################################################
# 3. Re-preform dimensionality reduction and clustering
################################################################################

# (1). Dimensionality reduction and clustering
B16OT1 <- NormalizeData(B16OT1)
B16OT1 <- FindVariableFeatures(B16OT1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(B16OT1)
B16OT1 <- ScaleData(B16OT1, features = all.genes)
B16OT1 <- RunPCA(B16OT1, features = VariableFeatures(object = B16OT1))
ElbowPlot(B16OT1)
B16OT1 <- FindNeighbors(B16OT1, dims = 1:20)
B16OT1 <- FindClusters(B16OT1, resolution = 0.4)
B16OT1 <- RunUMAP(B16OT1, dims = 1:20)

# Find Marker genes
B16OT1.filter.markers <- FindAllMarkers(B16OT1, only.pos = TRUE)
write.csv(B16OT1.filter.markers,file = '~/space/genes_oe/3_1.B16OT1.filter.markers.csv')

# Dimplots
pdf(file="~/space/genes_oe/3_1_dimplot.B16OT1.pdf",width = 5,height = 4.5)
DimPlot(B16OT1, reduction = "umap",label = T,cols = mypal_1)
dev.off()


# (2). UMAP split by sample
unique(B16OT1$Sample_Name)
B16OT1$Group <- B16OT1$Sample_Name
B16OT1$Group[B16OT1$Group %in% c("EV","Luc")]="Control"
unique(B16OT1$Group)

# Set levels
all_groups <- unique(B16OT1$Group)
other_groups <- sort(setdiff(all_groups, "Control"))
new_levels <- c("Control", other_groups)
B16OT1$Group <- factor(B16OT1$Group, levels = new_levels)

# Set color for samples
colors_sample
names(colors_sample)[39] <- "Control"
Idents(B16OT1) <- B16OT1$Group

# Plot
pdf(file="~/space/genes_oe/3_2_B16OT1.dimplot.samplesplit.pdf",width=7.5,height=8)
DimPlot(B16OT1, reduction = "umap",label = F,split.by = "Group",cols = colors_sample,ncol = 5,raster=F,pt.size = 0.02)+
  theme(
    plot.title = element_blank(),    
    strip.text = element_blank(),    
    strip.background = element_blank() 
  )
dev.off()


################################################################################
# 4. Pseudobulk PCA of each condition
################################################################################

# Calculate pseudobulk values
pb <- AggregateExpression(B16OT1, 
                          group.by = c("Group"), 
                          return.seurat = F,
                          normalization.method = NULL)
counts_pb <- pb$RNA
counts_pb <- as.data.frame(counts_pb)
meta_pb <- data.frame(group=colnames(counts_pb))
row.names(meta_pb) <- meta_pb$group
head(meta_pb)

# Construct DESeq2 DDS
dds <- DESeqDataSetFromMatrix(countData = counts_pb,
                              colData = meta_pb,
                              design = ~ group)

# Filter genes with non expression in all samples
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep, ]

# VST normalization for PCA
vst <- vst(dds, blind = TRUE)
vst_mat <- assay(vst)

# PCA + ggplot
pca <- prcomp(t(vst_mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>% 
  rownames_to_column("sample") %>% 
  left_join(meta_pb %>% rownames_to_column("sample"), by = "sample")

pca_df$group <- factor(pca_df$group,levels=c("Control",pca_df$group[-which(pca_df$group=="Control")]))

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4,shape=16) +
  scale_color_manual(values = colors_sample)+
  theme_bw() +
  theme(aspect.ratio = 1)+
  ggtitle("PCA Plot - Pseudobulk") +
  xlab(paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))
ggsave(filename = "~/space/genes_oe/4_1_B16OT1_pesudobulk.pca.pdf",pca_plot,width=5,height=5)


################################################################################
# 5. Calculate DEGs between the condition of each Gene Over-Expression and the Control group
################################################################################

# (1). Set parameters
B16OT1$Group_Condition <- as.character(B16OT1$Sample_Name)
B16OT1$Group_Condition[B16OT1$Group_Condition %in% c("EV","Luc")] <- "Control"

all_samples <- unique(B16OT1$Group_Condition)
target_samples <- setdiff(all_samples, c("Control"))
target_samples

# (2). Perform DEG analysis
DEGs_B16OT1 <- list()
for (sample in target_samples) {
  message(paste("Calculate:", sample, "vs Control..."))
  
  # ident.2 is Ctrl
  markers <- FindMarkers(
    B16OT1,
    ident.1 = sample,
    ident.2 = "Control",
    group.by = "Group_Condition",
    test.use = "wilcox" 
  )
  markers$gene <- row.names(markers)
  DEGs_B16OT1[[sample]] <- markers
}
names(DEGs_B16OT1)
head(DEGs_B16OT1[[1]])

for (sample in target_samples) {
  DEGs_B16OT1[[sample]]$neg10padj <- -log10( DEGs_B16OT1[[sample]]$p_val_adj)
  DEGs_B16OT1[[sample]]$neg10padj[is.infinite(DEGs_B16OT1[[sample]]$neg10padj)]=300
  
}
head(DEGs_B16OT1[[1]])

# Merge results into a big table
for (sample in names(DEGs_B16OT1)) {
  DEGs_B16OT1[[sample]]$sample <- paste0("",sample)
}

DEGs_B16OT1 <- do.call(rbind,DEGs_B16OT1)
head(DEGs_B16OT1)
write.csv(DEGs_B16OT1,file="~/space/genes_oe/5_1_degs_B16OT1.csv")

# Filter of significant DEGs
DEGs_B16OT1_sig <- DEGs_B16OT1[DEGs_B16OT1$p_val_adj<0.05,]
DEGs_B16OT1_sig <- DEGs_B16OT1_sig %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))
write.csv(DEGs_B16OT1_sig,file="~/space/genes_oe/5_1_degs_B16OT1_sig.csv")


################################################################################
# 6. hdWGCNA to identify gene co-expression modules
################################################################################

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

# Set the cowplot theme for ggplot
theme_set(theme_cowplot())
set.seed(42)

# Enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- SetupForWGCNA(
  B16OT1,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "B16OT1" 
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Group"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Group' # set the Idents of the metacell seurat object
)

# Normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' 
)

# Plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# Assemble with patchwork
pdf(file="~/space/genes_oe/6_1_hdWGCNA_PlotSoftPowers.pdf",width=8,height=6)
wrap_plots(plot_list, ncol=2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# Construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'B16OT1'
)

pdf(file="~/space/genes_oe/6_2_hdWGCNA_Dendrogram.pdf",width=8,height=4)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

# Calculate Module Eigengenes
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)

# Get Module Eigengenes
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# Compute Eigengenes-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "B16OT1-M"
)

# Visualize hub genes with the top kMEs
p <- PlotKMEs(seurat_obj, ncol=3)
pdf(file="~/space/genes_oe/6_3_hdWGCNA_KMEs.pdf",width=8,height=4)
print(p)
dev.off()

# Get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 20)
head(hub_df)

# Save hdWGCNA results
saveRDS(seurat_obj, file='~/space/genes_oe/6_3_hdWGCNA_onlyT.rds')



# Compute AUCell scoring for the top 20 hub genes by kME for each module
library(AUCell)
hub_list <- split(
  hub_df$gene_name,
  hub_df$module
)

sparse_expr <- GetAssayData(B16OT1, layer = "data")
# Run AUCell_buildRankings
cells_rankings_B16OT1 <- AUCell_buildRankings(
  sparse_expr,
  plotStats = FALSE,
  nCores = 8
)

# Run AUCell based on hub genes
cells_AUC_hubgene <- AUCell_calcAUC(
  hub_list, 
  cells_rankings_B16OT1,
  nCores = 8)

# Save AUCell results
saveRDS(cells_AUC_hubgene,"~/space/genes_oe/6_4_cells_AUC_hubgene.rds")


# Add AUCell results to seurat metadata
cells_AUC_hubgene <- as.data.frame(getAUC(cells_AUC_hubgene))
auc_matrix_t <- t(cells_AUC_hubgene)
all.equal(rownames(auc_matrix_t), rownames(B16OT1@meta.data))
B16OT1@meta.data <- cbind(B16OT1@meta.data, auc_matrix_t)
head(B16OT1@meta.data)


################################################################################
# 7. Visualize of gene co-expression modules results
################################################################################

# (1). Heatmap visualization of CM hub gene AUCell score differences between each gene over-expression condition and the control group
groups <- unique(B16OT1$Group)
groups <- setdiff(groups, "Control")
groups
pathways <- paste0("B16OT1-M",1:6)

# Create a result matrix
diff_matrix <- matrix(NA, 
                      nrow = length(groups), 
                      ncol = length(pathways),
                      dimnames = list(groups, pathways))

# Calculate the average difference between each Group and Ctrl (mean(Group) - mean(Ctrl))
ctrl_cells <- B16OT1$Group == "Control"

for (i in seq_along(groups)) {
  grp <- groups[i]
  grp_cells <- B16OT1$Group == grp
  
  for (j in seq_along(pathways)) {
    pw <- pathways[j]
    
    mean_grp <- mean(B16OT1@meta.data[[pw]][grp_cells], na.rm = TRUE)
    mean_ctrl <- mean(B16OT1@meta.data[[pw]][ctrl_cells], na.rm = TRUE)
    
    diff_matrix[i, j] <- mean_grp - mean_ctrl
  }
}
diff_matrix <- as.data.frame(diff_matrix)
diff_matrix$Group <- rownames(diff_matrix)
head(diff_matrix)

# Visualization
plot_matrix <- diff_matrix[, !colnames(diff_matrix) %in% "Group"]
rownames(plot_matrix) <- diff_matrix$Group

mat <- as.matrix(plot_matrix)
max_abs <- max(abs(mat), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
write.csv(plot_matrix,file="~/space/genes_oe/B16OT1/7_1_hdwgcna_AUCelldif.csv")

pdf(file="~/space/genes_oe/7_1_hdwgcna_AUCelldif_heatmap.pdf",width = 6,heigh=8)
pheatmap(plot_matrix,
         scale = "none",          
         cluster_rows = TRUE,     
         cluster_cols = TRUE,     
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Module Activity Dif vs Ctrl",
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = 45,           
         breaks = breaks,
         cellwidth = 10,
         cellheight = 10,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         border_color = "grey90",
         legend_breaks = seq(-max_abs, max_abs, length.out = 5))
dev.off()



# (2). Barplot visualization of AUCell score differences of CM5 between each gene over-expression condition and the control group
# Extract the M5 column and sort
df_plot <- data.frame(
  Gene = rownames(plot_matrix),
  Value = plot_matrix[, "B16OT1-M5"]
) %>%
  arrange(Value) 

df_plot$Gene <- factor(df_plot$Gene, levels = df_plot$Gene)

# Set color parameter
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# Visualization
p <- ggplot(df_plot, aes(x = Gene, y = Value, fill = Value)) +
  geom_col() + 
  coord_flip() +  
  scale_fill_gradientn(colors = my_palette, 
                       limits = c(-max(abs(df_plot$Value)), max(abs(df_plot$Value)))) +
  theme_minimal() +
  labs(title = "B16OT1-M5",
       x = "Samples",
       y = "Difference (M5 vs Ctrl)",
       fill = "Delta") +
  theme(
    axis.text.y = element_text(color = "black"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black")
  )

# Save the plot
ggsave("~/space/genes_oe/7_2_hdwgcna_AUCelldif_M5bar.pdf", 
       plot = p, width = 5, height = 8)



# (3). Featureplot visualization of AUCell scores of each CM in cells

# Set colors
features <- c("B16OT1-M1", "B16OT1-M2", "B16OT1-M3", "B16OT1-M4", "B16OT1-M5", "B16OT1-M6")
target_cols <- c("#61C0BA", "#EFEA3C", "#18499E", "#9E2C29", "#E7211A", "#6AB82D")

# Set a list for plots
plot_lists <- list()

for (i in 1:length(features)) {
  feature_data <- FetchData(B16OT1, vars = features[i])[[1]]
  min_val <- min(feature_data, na.rm = TRUE)
  max_val <- max(feature_data, na.rm = TRUE)
  mid_val <- (min_val + max_val) / 2
  p <- FeaturePlot(B16OT1, 
                   features = features[i], 
                   pt.size = 0.02, 
                   raster = F,
                   combine = FALSE,order = T)[[1]] +
    scale_color_gradient2(
      low = "grey75",
      mid = "grey95",
      high = target_cols[i],
      midpoint = mid_val,
      labels = c("-", "+"), 
      breaks = c(min_val, max_val),
      guide = guide_colorbar(
        ticks = FALSE,
        barwidth = 0.5,
        barheight = 4
      )
    )+
    theme_void() +
    theme(
      plot.title = element_blank(),
      legend.position = "right",
      aspect.ratio = 1,
      legend.title = element_blank() 
    )
  
  plot_lists[[i]] <- p
}

combined_plot <- wrap_plots(plot_lists, ncol = 2)
pdf(file="~/space/genes_oe/7_3_hdWGCNA_KMEs_featureplot_AUCell.pdf",width=4,heigh=6)
combined_plot
dev.off()


################################################################################
# 8. Calculate AUCell scores of STAR-T signature genes 
################################################################################

# (1). Load STAR-T signature genes and calculate AUCell scores
start_signature <- read.csv("~/space/analysis/domain/tcell/start/1_5_df_deg_cancer_target_up_final.csv") # STAR-T signature genes was generated in Step12, line 125
start_signature <- start_signature$Var1
start_signature

# Convert to mouse gene symbols
library("babelgene")
orthologs(genes = start_signature, species = "mouse")
start_signature <- orthologs(genes = start_signature, species = "mouse")$symbol
start_signature

# Calculate AUCell scores
cells_AUC_start_signature_B16OT1 <- AUCell_calcAUC(
  list("START_signature" = start_signature), 
  cells_rankings_B16OT1,
  nCores = 8)

# Save AUCell results
saveRDS(cells_AUC_start_signature_B16OT1,"~/space/genes_oe/8_1_cells_AUC_start_signature_B16OT1.rds")


# Add AUCell results to seurat metadata
cells_AUC_start_signature_B16OT1 <- as.data.frame(getAUC(cells_AUC_start_signature_B16OT1))
auc_matrix_t <- t(cells_AUC_start_signature_B16OT1)
all.equal(rownames(auc_matrix_t), rownames(B16OT1@meta.data))
B16OT1@meta.data <- cbind(B16OT1@meta.data, auc_matrix_t)
head(B16OT1@meta.data)


# (2). Calculate AUCell score differences between each gene over-expression groups and the control group
groups <- unique(B16OT1$Group)
groups <- setdiff(groups, "Control")
groups
pathways <- "START_signature"

diff_matrix <- matrix(NA, 
                      nrow = length(groups), 
                      ncol = length(pathways),
                      dimnames = list(groups, pathways))
ctrl_cells <- B16OT1$Group == "Control"

for (i in seq_along(groups)) {
  grp <- groups[i]
  grp_cells <- B16OT1$Group == grp
  
  for (j in seq_along(pathways)) {
    pw <- pathways[j]
    
    mean_grp <- mean(B16OT1@meta.data[[pw]][grp_cells], na.rm = TRUE)
    mean_ctrl <- mean(B16OT1@meta.data[[pw]][ctrl_cells], na.rm = TRUE)
    
    diff_matrix[i, j] <- mean_grp - mean_ctrl
  }
}

# Convert into data.frame
diff_matrix <- as.data.frame(diff_matrix)
diff_matrix$Group <- rownames(diff_matrix)

plot_matrix <- diff_matrix[, !colnames(diff_matrix) %in% "Group",drop=F]
rownames(plot_matrix) <- diff_matrix$Group

# (3). Barplot to visualize the AUCell score differences between each gene over-expression groups and the control group
df_plot <- data.frame(
  Gene = rownames(plot_matrix),
  Value = plot_matrix[, "START_signature"]
) %>%
  arrange(Value) 

df_plot$Gene <- factor(df_plot$Gene, levels = df_plot$Gene)
write.csv(df_plot,file="~/space/genes_oe/8_2_AUCelldif_start_signature.csv")

my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# Draw plots
p <- ggplot(df_plot, aes(x = Gene, y = Value, fill = Value)) +
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colors = c("#a1d8b1", "white", "#e69138"), 
                       limits = c(-max(abs(df_plot$Value)), max(abs(df_plot$Value)))) +
  theme_minimal() +
  labs(title = "START-Signature",
       x = "Samples",
       y = "Difference vs Ctrl",
       fill = "Delta") +
  theme(
    axis.text.y = element_text(color = "black"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black")
  )

ggsave("~/space/genes_oe/8_2_AUCelldif_start_signature.pdf", 
       plot = p, width = 5, height = 8)


################################################################################
# 9. Visualize of DEGs using volcano plots
################################################################################

# Prepare data
degs_target <- DEGs_B16OT1_sig[DEGs_B16OT1_sig$sample %in% c("Adgre5","Ccr5","Lipa","Tbx21","Zfp683"),]
plot_data <- degs_target
plot_data$sample <- factor(plot_data$sample, levels = names(colors_sample))


plot_data <- plot_data %>%
  mutate(
    is_sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
      TRUE ~ "Non-sig"
    )
  )

# Set colors
color_up <- "#D6604D"   
color_down <- "#4393C3"  
color_nonsig <- "lightgrey" 

# Plots
sample_list <- c("Adgre5", "Ccr5", "Lipa", "Tbx21", "Zfp683")
plot_list <- list()
for (s in sample_list) {
  
  # Extract data
  df_s <- plot_data %>% filter(sample == s)
  max_abs_fc <- max(abs(df_s$avg_log2FC), na.rm = TRUE)
  max_x <- ceiling(max_abs_fc)
  # Plot
  p_s <- ggplot(df_s, aes(x = avg_log2FC, y = neg10padj)) +
    
    # Signaficant level
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed",    
               color = "black",       
               linewidth = 0.5) +
    geom_vline(xintercept = 0, 
               linetype = "dashed",    
               color = "black",       
               linewidth = 0.5) +
    # First layer: unobtrusive gray dots
    geom_point(
      data = df_s %>% filter(is_sig == "Non-sig"),
      aes(x = avg_log2FC, y = neg10padj),
      color = color_nonsig, alpha = 0.7, size = 2,shape=16
    ) +
    
    # Second layer: significantly down-regulated blue dots (avg_log2FC < 0)
    geom_point(
      data = df_s %>% filter(is_sig == "Down"),
      aes(x = avg_log2FC, y = neg10padj),
      color = color_down, alpha = 0.7, size = 2,shape=16
    ) +
    
    # Third layer: significantly up-regulated blue dots (avg_log2FC > 0)
    geom_point(
      data = df_s %>% filter(is_sig == "Up"),
      aes(x = avg_log2FC, y = neg10padj),
      color = color_up, alpha = 0.7, size = 2,shape=16
    ) +
    xlim(-max_x, max_x) +
    theme_bw() +
    labs(
      title = paste("", s),
      x = "log2_FC",
      y = "-log10(p-adjust)"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      aspect.ratio = 1,
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      axis.text = element_text(color = "black")
    ) +
    
    # Add significant label
    geom_text_repel(
      data = df_s %>% filter(gene %in% c("Gzmb","Gzmc") & is_sig %in% c("Up", "Down")),
      aes(label = gene),
      min.segment.length = 0,
      segment.color = "black",
      segment.size = 0.5,
      size = 5, 
      color = "black"
    )
  
  plot_list[[s]] <- p_s
}

combined_p <- wrap_plots(plot_list, ncol = 3)

# Export
ggsave("~/space/genes_oe/9_1_B16OT1.volcanoplot.pdf", 
       plot = combined_p, width = 10, height = 6)

