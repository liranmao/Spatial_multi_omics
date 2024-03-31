rm(list=ls())
getwd()
setwd("/mnt/d/Users/Liran/38nanobody_deep/Spnanob38_deep_RNA")

# install.packages('Rfast2')
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(ArchR)
library(Rfast2)

source('./scripts/SpatialPlot_new.R')

sampleNames <- 'Spananob38_deep_RNA'
output_folder <- paste0(sampleNames, '_plot')

for (i in c(output_folder, paste0(output_folder,'/markers_plot/' ))){
  folder_path <- paste0("./",i)
  # Create the new folder
  dir.create(folder_path)
}



gene_matrix <- Read10X(data.dir = "/mnt/d/Users/Liran/38nanobody_deep/Spnanob38_deep_RNA/raw/")

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = gene_matrix, assay = assay)

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

p <- VlnPlot(spatial.obj, features = "nFeature_Spatial", pt.size = 0.1, log = TRUE) + NoLegend()
p
median(spatial.obj$nFeature_Spatial)

pdf(paste0('./', output_folder, '/', sampleNames, "_nfrags.pdf"), width = 5, height = 5)
print(p)
dev.off()


p <- SpatialPlot_new(spatial.obj, features = "nFeature_Spatial",  pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

pdf(paste0('./', output_folder, '/', sampleNames, "_nFrags_spatial.pdf"), width = 5, height = 5)
print(p)
dev.off()


#SpatialDimPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4)
#Idents(spatial.obj) <- 'Clusters'
#SpatialDimPlot(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = c('C1', 'C2')), facet.highlight = TRUE, ncol = 2, pt.size.factor = 3)


spatial.obj <- SCTransform(spatial.obj, assay = "Spatial", verbose = FALSE)

spatial.obj <- RunPCA(spatial.obj, assay = "SCT", verbose = FALSE)
spatial.obj <- FindNeighbors(spatial.obj, reduction = "pca", dims = 1:30)
spatial.obj <- FindClusters(spatial.obj, verbose = FALSE, resolution = 0.8)
spatial.obj <- RunUMAP(spatial.obj, reduction = "pca", dims = 1:30)

# change the color
n_clusters <- length(unique(spatial.obj$seurat_clusters))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0( seq_len(n_clusters)-1)
cols_swith <- cols
cols_swith[1] <- cols[6]
cols_swith[4] <- cols[1]
cols_swith[6] <- cols[4]

cols_swith[2] <- "#8A9FD1" 
cols_swith[7] <- "#272E6A" 

p1 <- DimPlot(spatial.obj, reduction = "umap", label = TRUE, cols = cols_swith)
pdf(paste0('./', output_folder, '/', sampleNames, "_clusters_umap.pdf"), width = 5, height = 4)
print(p1)
dev.off()

p2 <- SpatialDimPlot(spatial.obj, label = FALSE, label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols_swith) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p2$layers[[1]]$aes_params <- c(p2$layers[[1]]$aes_params, shape=22)

pdf(paste0('./', output_folder, '/', sampleNames, "_clusters_spatial.pdf"), width = 5, height = 5)
print(p2)
dev.off()


de_markers <- FindMarkers(spatial.obj, ident.1 = 7)
SpatialFeaturePlot(object = spatial.obj, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))

# search for features exhibiting spatial patterning in the absence of pre-annotation. 
spatial.obj <- FindSpatiallyVariableFeatures(spatial.obj, assay = "SCT", features = VariableFeatures(spatial.obj)[1:1000],
                                             selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(spatial.obj, selection.method = "moransi"), 10)


p3 <- SpatialFeaturePlot(spatial.obj, features = top.features, ncol = 1, alpha = c(0.1, 1),pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3

plot_features <- function(feature){
  p3 <- SpatialFeaturePlot(spatial.obj, features = feature, ncol = 1, alpha = c(0.1, 1),pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
  print(p3)
}

for (marker_gene in top.features){
  pdf(paste0('./', output_folder, '/markers_plot/', marker_gene, ".pdf"), width = 5, height = 5)
  plot_features(marker_gene)
  dev.off()
}

ggList <- lapply(top.features, plot_features)
ggList[[1]]

lapply(seq(length(ggList)),
       function(x)ggsave(filename=paste0('./', output_folder,'/markers_plot/', top.features[x],".png"), plot=ggList[[x]]))

saveRDS(spatial.obj, './seurat.rna.rds')


##################### Ignore the code below

markers_cluster <- as.data.frame(markerList_pos$C7)
markers_cluster$name[1:5]

features_spatial <- c('Gbx2', 'Prkcq', 'Smok4a', 'Cux2', 'Triqk', 
                      'Zfp536', 'Tenm4')
#features_spatial <- markerGenes
# feature <- features_spatial[1]
#feature <- 'Fgfr3'

# p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
#   theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
# p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
# p
# 
# png(filename = paste0('./markers_plot/', feature, '_upReg_spatial.png'), width = 1200, height = 1200, res = 300)
# p
# dev.off()


# plot list
# features_spatial <- markerGenes
#features_spatial <- as.data.frame(markerList_pos$C7)$name[11:13]

plot_features <- function(feature){
  p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

# source('SpatialPlot_new_set_limit.R')
# plot_features <- function(feature){
#   p <- SpatialPlot_new_set_limit(spatial.obj, features = feature, pt.size.factor = 4, image.alpha = 0, stroke = 0, limits=c(0.32,0.55)) +
#     theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
#   p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
#   p
# }

spatial.obj@assays$Spatial
saveRDS(spatial.obj, 'D:/Users/Liran/38nanobody/combine_38nanobody/seurat.rna.rds')
saveRDS(seurat.wnn.peak, 'D:/Users/Liran/38nanobody/combine_38nanobody/seurat.wnn.peak.new_with_all.rds')


spatial.obj
query
