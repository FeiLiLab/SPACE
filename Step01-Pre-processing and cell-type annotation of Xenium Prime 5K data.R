library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ggplot2)
library(Cairo)
library(arrow)
library(ggrepel)
library(BiocParallel)
library(future)
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
options(future.globals.maxSize = 300 * 1024^3)

################# 1. Construct a Seurat object from Xenium assay and perform QC
################# Bladder cancer data were used as an example.
#################
#################
#################
#################

################################################################################
# 1_1. Load 10x Genomics Xenium in-situ data
################################################################################

obj_BLCA <- LoadXenium("~/space/rawdata/BLCA/", fov = "fov")
obj_BLCA

################################################################################
# 1_2. Identify high-quality cells
################################################################################

obj_BLCA <- subset(obj_BLCA, subset = nCount_Xenium > 30 & nFeature_Xenium > 10)
obj_BLCA

################################################################################
# 1_3. Visualize the spatial distribution of cells
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/1_3_ImageDimPlot_sample.png",width=1900,height=900)
ImageDimPlot(obj_BLCA,size = 0.15,axes = TRUE,dark.background = FALSE) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
        plot.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))  
dev.off()

################################################################################
# 1_4. Visualize the expression of lineage markers using the ImageDimPlot function
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/1_4_ImageDimPlot_lineagemarker.png",width=1900,height=900)
ImageDimPlot(obj_BLCA,molecules = c("EPCAM","COL5A1","PECAM1","CD68"),size = 0.15,axes = TRUE,dark.background = FALSE) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
        plot.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))
dev.off()

################################################################################
# 1_5.  Visualize the expression of lineage markers using the ImageFeaturePlot function
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/1_5_ImageFeaturePlot_lineagemarker.png",width=1900,height=900)
options(repr.plot.width=16, repr.plot.height=10) 
ImageFeaturePlot(obj_BLCA, features = c("EPCAM","COL5A1","PECAM1","CD68"), 
                 max.cutoff = c(5,5,5,5), size = 0.75, cols = c("grey95", "red"),dark.background = FALSE)
dev.off()

################# 2. Annotate sample information
#################
#################
#################
#################
#################
#################

################################################################################
# 2_1. Extract coordinates for all the samples
################################################################################

library(dplyr)
library(purrr)
setwd("~/space/analysis/BLCA/coordinate/")
file_paths <- list.files(pattern = "BLCA_\\d+_coordinates\\.csv") %>%
  stringr::str_sort(numeric = TRUE)
process_file <- function(file_path) {
  sample_name <- gsub("_coordinates\\.csv", "", basename(file_path))
  coord_data <- tryCatch({
    raw_data <- read.csv(file_path, row.names = NULL)[-c(1,2), ]
    colnames(raw_data) <- c("X", "Y")
    data.frame(
      Sample = sample_name,
      X_max = max(as.numeric(raw_data$X)),
      X_min = min(as.numeric(raw_data$X)),
      Y_max = max(as.numeric(raw_data$Y)),
      Y_min = min(as.numeric(raw_data$Y))
    )
  }, error = function(e) {
    warning(paste("Error processing", file_path, ":", e$message))
    data.frame(Sample = sample_name,
               X_max = NA, X_min = NA,
               Y_max = NA, Y_min = NA)
  })
  return(coord_data)
}
result_df <- map_dfr(file_paths, process_file)
head(result_df)
dim(result_df)
class(result_df)
write.csv(result_df, "~/space/analysis/BLCA/2_1_BLCA_coordinates_summary.csv", row.names = FALSE)

################################################################################
# 2_2. Add coordinate information to meta.data
################################################################################

coords <- GetTissueCoordinates(obj_BLCA, which = "centroids")
head(coords)
obj_BLCA$coord_x<-coords$x
obj_BLCA$coord_y<-coords$y
head(obj_BLCA@meta.data)

################################################################################
# 2_3. Assign cells to samples
################################################################################

meta.data<-obj_BLCA@meta.data
sample<-result_df$Sample
cell_sample<-data.frame()
for (i in 1:length(sample)){
  x_max <-result_df[which(result_df$Sample==sample[i]),]$X_max
  x_min <-result_df[which(result_df$Sample==sample[i]),]$X_min
  y_max <-result_df[which(result_df$Sample==sample[i]),]$Y_max
  y_min <-result_df[which(result_df$Sample==sample[i]),]$Y_min
  
  cell_select<- rownames(meta.data[which(
    meta.data$coord_x >= x_min &
      meta.data$coord_x <= x_max &
      meta.data$coord_y >= y_min &
      meta.data$coord_y <= y_max),])
  cell_sample1<-data.frame(cell=cell_select,
                           sample=rep(sample[i]),length(cell_select))
  cell_sample<-rbind(cell_sample,cell_sample1)
}
dim(cell_sample)
length(unique(cell_sample$cell))
write.csv(cell_sample,"~/space/analysis/BLCA/2_3_BLCA_cell_sample.csv")
head(cell_sample)

################################################################################
# 2_4. Check if there are overlapping cells between samples
################################################################################

cell_sample_dup <- cell_sample %>% group_by(cell) %>%filter(n() > 1)
table(cell_sample_dup$sample)


################################################################################
# 2_5. Extract cells (this step will filter out unassigned cells)
################################################################################

obj_BLCA <- subset(obj_BLCA, cells = cell_sample$cell)
obj_BLCA@meta.data[1,]

row.names(cell_sample)<-cell_sample$cell
cell_sample<-cell_sample[row.names(obj_BLCA@meta.data),]
obj_BLCA$sample<-cell_sample$sample
obj_BLCA$sample<-factor(obj_BLCA$sample,levels=names(table(obj_BLCA$sample)))

CairoPNG("~/space/analysis/figure/BLCA/2_5_ImageDimPlot_sampleinfo.png",width=1900,height=900)
ImageDimPlot(obj_BLCA,size = 0.15,cols=mypal_1,group.by="sample",axes = TRUE,dark.background = FALSE) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
        plot.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))  
dev.off()


################# 3. Perform normalization, dimensionality reduction and clustering
#################
#################
#################
#################
#################
#################

options(future.globals.maxSize = 300 * 1024^3)
obj_BLCA <- NormalizeData(obj_BLCA, normalization.method = "LogNormalize", scale.factor = 100)
obj_BLCA <- ScaleData(obj_BLCA)
obj_BLCA <- RunPCA(obj_BLCA,npcs = 30, features = rownames(obj_BLCA))
obj_BLCA <- FindNeighbors(obj_BLCA, dims = 1:30,k.param = 20,features = rownames(obj_BLCA))  
obj_BLCA <- RunUMAP(obj_BLCA, dims = 1:30)
obj_BLCA <- FindClusters(obj_BLCA, resolution = 0.5)
Idents(obj_BLCA)<-obj_BLCA$Xenium_snn_res.0.5
Idents(obj_BLCA)<-factor(Idents(obj_BLCA),levels=c(0:length(unique(Idents(obj_BLCA)))))
table(Idents(obj_BLCA))
head(obj_BLCA@meta.data)

################# 4. Perform visualization
#################
#################
#################
#################
#################
#################

################################################################################
# 4_1. Visualize cell clustering
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/4_1_DimPlot_cluster.png",width=1000,height=900)
DimPlot(obj_BLCA,cols=mypal_3,raster=FALSE)
dev.off()

CairoPNG("~/space/analysis/figure/BLCA/4_1_DimPlot_cluster_sample.png",width=1000,height=900)
DimPlot(obj_BLCA,group.by="sample",raster=FALSE)
dev.off()

CairoPNG("~/space/analysis/figure/BLCA/4_1_DimPlot_inte_reos0.5.png",width=1000,height=900)
DimPlot(obj_BLCA,cols=mypal_3,raster=FALSE)
dev.off()

################################################################################
# 4_2. Visualize QC information-FeaturePlot 
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/4_2_FeaturePlot_QC.png",width=1800,height=500)
FeaturePlot(obj_BLCA, features = c("nFeature_Xenium", "nCount_Xenium"),raster=TRUE,min.cutoff =0,ncol=4,reduction = "umap",cols= c("gray90", "red"))
dev.off()

################################################################################
# 4_3. Visualize QC information-Vlnplot
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/4_3_VlnPlot_QC.png",width=800,height=500)
VlnPlot(obj_BLCA, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()

################################################################################
# 4_4. Visualize the expression of marker gene
################################################################################

markerxe_homo<-intersect(c("EPCAM","CDH1", # Malignant
                           "CD3E","CD8A","CD4",#T_cell
                           "MS4A1","CD79A",#B_cell
                           "XBP1","MZB1",#Plasma
                           "FCGR2A","CD14",#Myeloid
                           "NOTCH3","PDGFRB", #Pericyte
                           "COL5A1","COL5A2",#Fibroblast
                           "CSRP1","MYLK",#Smooth_muscle
                           "PECAM1","CD34",#Endothelial
                           "MPZ","SCN7A"#Neural-like
),row.names(obj_BLCA))

CairoPNG("~/space/analysis/figure/BLCA/4_4_FeaturePlot_Marker.png",width=2000,height=3000)
FeaturePlot(obj_BLCA, features = markerxe_homo,raster=TRUE,min.cutoff =0,ncol=4,reduction = "umap",cols= c("gray90", "red"))
dev.off()

################################################################################
# 4_5. Visualize cell clustering by sample
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/4_5_DimPlot_sample.png",width=1200,height=900)
DimPlot(obj_BLCA,group.by="sample",raster=FALSE)
dev.off()

################################################################################
# 4_6. Identify high-quality cell clusters (this step will filter out the low-quality cell clusters with low cell number and low gene number detected)
################################################################################

cluster_stats <- as.data.frame(obj_BLCA@meta.data %>%
                                 group_by(seurat_clusters) %>%
                                 summarise(
                                   mean_nFeature = mean(nFeature_Xenium),
                                   mean_nCount = mean(nCount_Xenium),
                                   n_cells = n()))
cluster_stats_retained<-cluster_stats[(cluster_stats$n_cells>2&cluster_stats$mean_nFeature>100),]

################################################################################
# 4_7. Retain high-quality cell clusters for downstream analysis
################################################################################

obj_BLCA<-subset(obj_BLCA,seurat_clusters %in% as.character(cluster_stats_retained$seurat_clusters))
table(Idents(obj_BLCA))

################# 5. Perform cell annotation
#################
#################
#################
#################
#################
#################

plan("multisession",workers=2)

################################################################################
# 5_1. Identify marker genes for each cluster
################################################################################

markers_BLCA <- FindAllMarkers(obj_BLCA,only.pos = TRUE)
write.csv(markers_BLCA,"~/space/analysis/figure/BLCA/5_1_Markers_BLCA_reso0.5.csv")
table(Idents(obj_BLCA))
table(obj_BLCA$Xenium_snn_res.0.5)

################################################################################
# 5_2. Annotate cell clusters according to lineage markers (XXX represents cluster ID)
################################################################################

obj_BLCA <- RenameIdents(obj_BLCA,'XXX' = 'Malignant','XXX' = 'Fibroblast','XXX' = 'Smooth_muscle','XXX' = 'Endothelial','XXX' = 'Pericyte',
                         'XXX' = 'T_cell','XXX' = 'B_cell','XXX' = 'Plasma','XXX' = 'Myeloid','XXX' = 'Neural')
table(Idents(obj_BLCA))
obj_BLCA$celltype<-Idents(obj_BLCA)
table(obj_BLCA$celltype)

################################################################################
# 5_3. Perform visualization
################################################################################

CairoPNG("~/space/analysis/figure/BLCA/5_3_DimPlot_celltype.png",width=1000,height=900)
DimPlot(obj_BLCA,cols=mypal_3,raster=FALSE)
dev.off()

CairoPNG("~/space/analysis/figure/BLCA/5_3_ImageDimPlot_celltype.png",width=1900,height=900)
ImageDimPlot(obj_BLCA,size = 0.15,cols=mypal_3,axes = TRUE,dark.background = FALSE) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
        plot.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))  
dev.off()

CairoPNG("~/space/analysis/figure/BLCA/5_3_ImageDimPlot_celltype_sample.png",width=1900,height=900)
ImageDimPlot(obj_BLCA,size = 0.15,group.by="sample",cols=mypal_1,axes = TRUE,dark.background = FALSE) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),
        plot.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))  
dev.off()

pdf("~/space/analysis/figure/BLCA/5_3_cell.prop.sample.pdf",width=10,height=5)
cell.prop.sample<-as.data.frame(prop.table(table(Idents(obj_BLCA),obj_BLCA$sample)))
colnames(cell.prop.sample)<-c("cluster","sample","proportion")
ggplot(cell.prop.sample,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=mypal_3)
dev.off()

################################################################################
# 5_4. Save data
################################################################################

# Add cancer type suffix to cell IDs
obj_BLCA <- RenameCells(obj_BLCA, new.names = paste0(colnames(obj_BLCA), "_BLCA"))

metadata<-obj_BLCA@meta.data
saveRDS(metadata,"~/space/metadata/obj_BLCA_metadata.rds")
saveRDS(obj_BLCA,"~/space/rds/obj_BLCA.rds")
