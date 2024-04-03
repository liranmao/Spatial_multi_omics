rm(list=ls())

library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

source('./scripts/SpatialPlot_new.R')

sampleNames <- 'Spnbprot68'

# read ADT
adt <- Read10X(data.dir = "./umi_count/", gene.column=1)
adt <- adt[1:(nrow(adt) - 1), , drop = FALSE]
adt

###
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

#object <- CreateSeuratObject(counts = merged_matrix, assay = assay)
object <- CreateSeuratObject(counts = adt, assay = assay)


image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

p <- VlnPlot(spatial.obj, features = "nFeature_Spatial", pt.size = 0.1, log = TRUE) + NoLegend()
p
median(spatial.obj$nFeature_Spatial)
# png(filename = 'nFrags.png', width = 1200, height = 1200, res = 300)
# p
# dev.off()

UMI_count <- rowSums(spatial.obj[,2:(ncol(spatial.obj)-1)])
#Count Proteins per pixel
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
protein_count <- rowSums(data_filtered_binary)

# for adt normalization
spatial.obj <- NormalizeData(spatial.obj, normalization.method = "CLR", margin = 2)
spatial.obj <- RunPCA(spatial.obj, assay = "CLR", verbose = FALSE)

## Plot ADT
SpatialFeaturePlot(object = spatial.obj, features = c('CD3-GTATGTCCGCTCGAT'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q99", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
SpatialFeaturePlot(object = spatial.obj, features = c('CD4-AACAAGACCCTTGAG'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q5", max.cutoff = "q99", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
SpatialFeaturePlot(object = spatial.obj, features = c('CD140a-GTCATTGCGGTCCTA'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q15", max.cutoff = "q90", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
SpatialFeaturePlot(object = spatial.obj, features = c('CD133-CTAGACCCTTCCCTT'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q15", max.cutoff = "q97", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
SpatialFeaturePlot(object = spatial.obj, features = c('B220-CCTACACCTCATAAT'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q5", max.cutoff = "q97", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
SpatialFeaturePlot(object = spatial.obj, features = c('CD90.1-AGTATGGGATGCAAT'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q5", max.cutoff = "q97", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
SpatialFeaturePlot(object = spatial.obj, features = c('CD90.2-CCGATCAGCCGTTTA'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q5", max.cutoff = "q90", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))

SpatialFeaturePlot(object = spatial.obj, features = c('CD34-GATTCCTTTACGAGC'), alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q5", max.cutoff = "q90", image.alpha = 0.5, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))

