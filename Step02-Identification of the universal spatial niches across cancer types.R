library(Seurat)
library(SeuratWrappers)
library(dplyr)
options(future.globals.maxSize = 300 * 1024^3)


################# 1. Perform spatial niche identification with Seurat BuildNicheAssay function
#################
#################
#################
#################
#################


################################################################################
# 1_1. Data preparation for spatial niche identification
################################################################################

# (1). Load meta data
cancer_list <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")

for(i in cancer_list){
  assign(paste0(i,"_metadata"), readRDS(paste0("~/space/metadata/obj_",i,"_metadata.rds"))) # Meta data of each cancer were generated in Step01, line 376
  }

for (i in cancer_list) {
  print(eval(parse(text = paste0("range(",i,"_metadata$coord_x)"))))
  print(eval(parse(text = paste0("range(",i,"_metadata$coord_y)"))))
}

# (2). Adjust cell coordinates to avoid spatial overlap between slices
# All cells from different cancer slices were merged into a single unified coordinate system to enable simultaneous inference of spatial niches (SNs) across cancers
# To prevent spatial overlap between slices, we applied a fixed offset of +12,000 to the x-coordinates and +24,000 to the y-coordinates of each cancer slice
# These offset values were chosen based on the typical range of x and y coordinates of cancer slices for the Xenium assay, where x and y coordinates generally fall between 0 to 12,000, and between 0 to 24,000, respectively

# Set the offset values of cell coordinates of each cancer 
x_offset_vec <- c(
  BLCA = 0, BRCA = 12000, CRC = 24000, ESCC = 36000, GBM = 48000, LIHC = 60000,
  NSCLC = 0, PAAD = 12000, PRAD = 24000, RCC = 36000, SKCM = 48000, STAD = 60000
)

y_offset_vec <- c(
  BLCA = 0, BRCA = 0, CRC = 0, ESCC = 0, GBM = 0, LIHC = 0,
  NSCLC = 24000, PAAD = 24000, PRAD = 24000, RCC = 24000, SKCM = 24000, STAD = 24000
)


# Apply coordinate adjustments for each cancer
for (current_cancer in cancer_list) {

  x_off <- x_offset_vec[current_cancer]
  y_off <- y_offset_vec[current_cancer]
  
  obj_name <- paste0(current_cancer, "_metadata")

  current_obj <- get(obj_name)
  current_obj$coord_x_adj = current_obj$coord_x
  current_obj$coord_y_adj = current_obj$coord_y
  if (x_off != 0) current_obj$coord_x_adj <- current_obj$coord_x_adj + x_off
  if (y_off != 0) current_obj$coord_y_adj <- current_obj$coord_y_adj + y_off
  assign(obj_name, current_obj)
}


# (3). Combine metadata of all cancer types
list_of_dfs <- list(BLCA_metadata,BRCA_metadata,CRC_metadata,ESCC_metadata,GBM_metadata,LIHC_metadata,
                    NSCLC_metadata,PAAD_metadata,PRAD_metadata,RCC_metadata,SKCM_metadata,STAD_metadata)

allcancer_meta_merged <- do.call("rbind", list_of_dfs)

################################################################################
# 1_2. Perform spatial niche identification with Seurat BuildNicheAssay function
################################################################################
# We used the BuildNicheAssay function in the Seurat package to identify spatial niches
# Notice: To ensure compatibility with our metadata structure (rather than a full Seurat object), we adapted the function's original source code

# (1). Set parameters
niches.k = 5
neighbors.k = 20

# (2). Initialize an empty cells × groups binary matrix
cells <- row.names(allcancer_meta_merged)
group.labels_all <- unlist(allcancer_meta_merged["celltype"][cells, ])
groups <- sort(unique(group.labels_all))
cell.type.mtx <- matrix(
  data = 0,
  nrow = length(cells),
  ncol = length(groups)
)
rownames(cell.type.mtx) <- cells
colnames(cell.type.mtx) <- groups

# (3). Fill the binary matrix
cells.idx <- seq_along(cells)
group.idx <- match(group.labels_all, groups)
cell.type.mtx[cbind(cells.idx, group.idx)] <- 1

# (4). Find neighbors based on cell position
coords <- allcancer_meta_merged[,c("coord_x_adj","coord_y_adj")]
colnames(coords)=c("x","y")
coords <- as.matrix(coords[ , c("x", "y")])
neighbors <- FindNeighbors(coords, k.param = neighbors.k, compute.SNN = FALSE)

# (5). Create an assay for niche identification
sum.mtx <- as.matrix(neighbors[["nn"]] %*% cell.type.mtx)
niche.assay <- CreateAssayObject(counts = t(sum.mtx))

counts <- GetAssayData(niche.assay, slot = "counts")
seurat_obj <- CreateSeuratObject(counts = counts,
                                 project = "niche_proj",
                                 assay = "niche",       
                                 min.cells = 0, 
                                 min.features = 0)

DefaultAssay(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, layer = "counts")

# (6). Perform k-means to identify spatial niches
results <- kmeans(
  x = t(seurat_obj@assays$niche@layers$scale.data),
  centers = niches.k,
  nstart = 30
)

# (7). Append spatial niche information to the metadata
seurat_obj@meta.data$niche <- results[["cluster"]]
seurat_obj=AddMetaData(seurat_obj,allcancer_meta_merged)

# (8). Append spatial niche information to the RDS of each cancer
cancer_list <- c("BLCA", "BRCA", "CRC", "ESCC","GBM","LIHC", "NSCLC", "PAAD", "PRAD", "RCC", "SKCM","STAD")
for(cancer in cancer_list){
  obj <- readRDS(paste0("~/space/rds/obj_", cancer, ".rds")) # RDS data files were generated in Step01, line 377
  obj$niche <- seurat_obj@meta.data$niche[match(row.names(obj@meta.data), row.names(seurat_obj@meta.data))]
  message("Processed: ", cancer)
  saveRDS(obj,file=paste0("~/space/rds/obj_", cancer, ".rds"))
}

# (9). Export final results
saveRDS(seurat_obj@meta.data,file="~/space/seurat_niche/metadata_merge.rds")

