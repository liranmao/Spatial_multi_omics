setwd('/mnt/d/Users/Liran/17nanobody/combine_17nanobody')
samplepath_list <- c('/mnt/d/Users/Liran/17nanobody/Spnanob17_H3K27me3/outs/Save-inTissue-Spananob17_BC6',
                     '/mnt/d/Users/Liran/17nanobody/Spnanob17_H3K27ac/Save-inTissue-Spananob17_BC7')
modalitylist <- c('BC6_h3k27me3', 'BC7_h3k27ac')
modalities <- c('BC6_h3k27me3', 'BC7_h3k27ac')
colnameslist<-c('Spananob17_BC6#','Spananob17_BC7#')
outputname <-'Spananob17_combine'

for (i in c('combine_plot')){
  folder_path <- paste0("./",i)
  # Create the new folder
  dir.create(folder_path)
  
}


getAvailableMatrices(projCUTA)
test <- getMatrixFromProject( ArchRProj = projCUTA, useMatrix = "TileMatrix",binarize = TRUE)
test
getMatrixFromProject( ArchRProj = projCUTA, useMatrix = "PeakMatrix")

####### change name
# if only creat peak object
for (i in 1:2){
  samplepath <- samplepath_list[i]
  modality <- modalitylist[i]
  colname <- colnameslist[i]
  projCUTA <- loadArchRProject(path = samplepath, force = FALSE, showLogo = TRUE)
  
  projCUTA <- addGroupCoverages(ArchRProj = projCUTA, groupBy = "Clusters")
  pathToMacs2 <- '/home/liranmao/.local/bin/macs2'
  projCUTA <- addReproduciblePeakSet(
    ArchRProj = projCUTA, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2,
    force = TRUE
  )
  
  projCUTA <- addPeakMatrix(projCUTA)
  PeakSet <- getPeakSet(projCUTA)
  SummarizedExperiment<- getMatrixFromProject( ArchRProj = projCUTA, useMatrix = "PeakMatrix")
  SummarizedExperiment_assay <- assay(SummarizedExperiment)
  
  
  # Combine seqnames and ranges into column names
  col_names <- paste(seqnames(PeakSet), ranges(PeakSet), sep = "-")
  # Set the column names of the sparse matrix
  rownames(SummarizedExperiment_assay) <- col_names
  
  colnames_sparse <- colnames(SummarizedExperiment_assay)
  new_colnames <- sub(colname, "", colnames_sparse)
  new_colnames <- sub("-1", "", new_colnames)
  colnames(SummarizedExperiment_assay) <-new_colnames
  if(i == 1){
    seurat.wnn.peak <- CreateSeuratObject(counts = SummarizedExperiment_assay, assay = paste0('peaks_',modality))
  } else {
    # create a new assay to store ADT information
    adt_assay <- CreateAssayObject(counts = SummarizedExperiment_assay)
    # add this assay to the previously created Seurat object
    seurat.wnn.peak[[paste0('peaks_',modality)]] <- adt_assay
  }
  
}



# Save the Seurat object to an RDS file
saveRDS(seurat.wnn.peak, file = "/mnt/d/Users/Liran/17nanobody/combine_17nanobody/seurat.wnn.peak.rds")



######################
####### the following is all in peaks 

for(x in modalities){
  DefaultAssay(seurat.wnn.peak) <- paste0('peaks_',x)
  seurat.wnn.peak <- RunTFIDF(seurat.wnn.peak) %>% FindTopFeatures() %>% RunSVD(reduction.name = paste0(x,'_lsi'))
}


seurat.wnn.peak <- FindMultiModalNeighbors(
  seurat.wnn.peak, reduction.list = list( 'BC6_h3k27me3_lsi', 'BC7_h3k27ac_lsi'), 
  dims.list = list(2:30,2:30), modality.weight.name = "histone.weight",k.nn = 10
)

seurat.wnn.peak <- RunUMAP(seurat.wnn.peak, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat.wnn.peak <- FindClusters(seurat.wnn.peak, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)

seurat.wnn.peak <- FindClusters(seurat.wnn.peak, graph.name = "wsnn", algorithm = 3, resolution = 0.6, verbose = FALSE)
seurat.wnn.peak$wsnn_res.0.6

for(x in modalities){
  assay <- paste0('peaks_',x,'.weight')
  p2 <- FeaturePlot(seurat.wnn.peak,assay,min.cutoff = 0.2,max.cutoff = 0.6) + scale_color_viridis_c()
  print(p2)
}

n_clusters <- length(unique(seurat.wnn.peak$wsnn_res.0.8))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0( seq_len(n_clusters)-1)
cols

p3 <- DimPlot(seurat.wnn.peak,label=TRUE,pt.size=0.5, cols = cols)

pdf(paste0("combine_plot/", outputname, "_clusters_umap_wsnn_res.0.8.pdf"), width = 5, height = 4)
print(p3)
dev.off()

head(Idents(seurat.wnn.peak), 5)




## add image
data.dir <- "/mnt/d/Users/Liran/17nanobody/Spnanob17_H3K27ac"
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[sub("-1", "", sub("Spananob18_BC5#", "", Cells(x = seurat.wnn.peak)))]
DefaultAssay(seurat.wnn.peak = image) <- 'peaks_BC7_h3k27ac'
seurat.wnn.peak[[slice]] <- image


# spatial plot weight cluster
for(x in modalities){
  assay <- paste0('peaks_',x,'.weight')
  # p3 <- SpatialPlot_new(seurat.wnn.peak, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
  #   theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  # p3$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  # print(p)
  
  p3 <- SpatialPlot(seurat.wnn.peak, label = FALSE, label.size = 3, features = assay,pt.size.factor = 4.5,  image.alpha = 0, stroke = 0, min.cutoff = 0.3, max.cutoff = 0.9)+ scale_color_viridis_c()
  p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
  print(p3)
}



# spatial plot cluster
p4 <- SpatialPlot(seurat.wnn.peak, label = FALSE, label.size = 3, group.by = 'wsnn_res.0.8', pt.size.factor = 4.5, cols = cols, image.alpha = 0, stroke = 0)
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
p4

pdf(paste0("combine_plot/", outputname, "_clusters_spatial_wsnn_res.0.8.pdf"), width = 5, height = 5)
print(p4)
dev.off()

# add gene score
for (modality in modalitylist){
  gene_score <- readRDS(paste0('Spananob17_', modality,'_gene_score.rds'))
  # seurat <- AddAssay(seurat, counts = your_ATAC_data, assay = "ATAC")
  # create a new assay to store ADT information
  adt_assay <- CreateAssayObject(counts = gene_score)
  # add this assay to the previously created Seurat object
  seurat.wnn.peak[[paste0('gscore_',modality)]] <- adt_assay
}

# subcluster of liver and heart
# take liver as an example first
liver <- seurat.wnn.peak[,seurat.wnn.peak@meta.data$wsnn_res.0.8 %in% c(5)]

seurat.wnn.peak@meta.data$

p4<-SpatialDimPlot(liver, crop = TRUE, label = TRUE,label.size = 3,  pt.size.factor = 4.5, cols = cols, image.alpha = 1, stroke = 0)
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
p4

p3 <- DimPlot(liver,label=TRUE,pt.size=0.5, cols = cols)
p3

cluster6.markers <- FindMarkers(seurat.wnn.peak,assay = 'peaks_BC7_h3k27ac', ident.1 = 6, logfc.threshold = 0.25, test.use = "wilcox")

liver

subcluster <- function(inpute_seraut){
  liver_new <- DietSeurat(inpute_seraut, dimreducs = NULL)
  VariableFeatures(liver_new) <- names(which(Matrix::rowSums(liver_new) > 1))
  
  for(x in modalities){
    DefaultAssay(liver_new) <- paste0('peaks_',x)
    liver_new <- RunTFIDF(liver_new) %>% FindTopFeatures() %>% RunSVD(reduction.name = paste0(x,'_lsi'))
  }
  
  
  liver_new <- FindMultiModalNeighbors(
    liver_new, reduction.list = list( 'BC6_h3k27me3_lsi', 'BC7_h3k27ac_lsi'), 
    dims.list = list(2:30,2:30), modality.weight.name = "histone.weight",k.nn = 10, knn.range = (dim(liver_new@meta.data)[1]-1) # number of cells
  )
  
  liver_new <- RunUMAP(liver_new, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  liver_new <- FindClusters(liver_new, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)
  
  # liver_new <- FindClusters(liver_new, graph.name = "wsnn", algorithm = 3, resolution = 0.6, verbose = FALSE)
  # seurat.wnn.peak$wsnn_res.0.6
  # 
  
  p3 <- DimPlot(liver_new,label=TRUE,pt.size=0.5, cols = cols)
  print(p3)
  
  p4<-SpatialDimPlot(liver_new, crop = TRUE, label = TRUE,label.size = 3,  pt.size.factor = 4.5, cols = cols, image.alpha = 1, stroke = 0)
  p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
  print(p4)
  return(liver_new)
}

liver_new2 <- liver_new
liver_new2
heart <- seurat.wnn.peak[,seurat.wnn.peak@meta.data$wsnn_res.0.8 %in% c(3)]

p4<-SpatialDimPlot(heart, crop = TRUE, label = TRUE,label.size = 3,  pt.size.factor = 4.5, cols = cols, image.alpha = 1, stroke = 0)
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
p4

p3 <- DimPlot(heart,label=TRUE,pt.size=0.5, cols = cols)
p3

heart_new <- subcluster(heart)


liver_more <- seurat.wnn.peak[,seurat.wnn.peak@meta.data$wsnn_res.0.8 %in% c(5, 9, 12)]

p4<-SpatialDimPlot(liver_more, crop = TRUE, label = TRUE,label.size = 3,  pt.size.factor = 4.5, cols = cols, image.alpha = 1, stroke = 0)
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
p4

p3 <- DimPlot(liver_more,label=TRUE,pt.size=0.5, cols = cols)
p3



# time






##########haven't done things below yet
# go for a certain cluster
# get features, mimic archr, use https://satijalab.org/seurat/reference/findallmarkers, https://satijalab.org/seurat/articles/pbmc3k_tutorial
# get all markers for all cluster
markersCSS_combine <- FindAllMarkers(seurat.wnn.peak, assay = 'gscore_BC6_h3k27me3', logfc.threshold = 0.25, test.use = "wilcox")
markersGAS_combine <- FindAllMarkers(seurat.wnn.peak, assay = 'gscore_BC7_h3k27ac', logfc.threshold = 0.25, test.use = "wilcox")

# venn plot for genes
install.packages("VennDiagram")
# Load the VennDiagram package
library(VennDiagram)

css_df <- data.frame(Gene = markersCSS_combine$gene, Sample = 'CSS')
gas_df <- data.frame(Gene = markersGAS_combine$gene, Sample = 'GAS')
venn_df <- rbind(css_df,gas_df)
write.csv(venn_df, "css_gas_gene_score_venn_df.csv", row.names = FALSE)


# Create a Venn diagram
venn.plot <- venn.diagram(
  x = list(list1 = markersCSS_combine$gene, list2 = markersGAS_combine$gene),
  category.names = c("markersCSS", "markersGAS"),
  filename = NULL
)

# Display the Venn diagram
grid.draw(venn.plot)

intersection_gene <- intersect(markersCSS_combine$gene, markersGAS_combine$gene)
print(intersection_gene)
length(intersection_gene)
'Nfe2'%in% intersection_gene
'Hand2' %in% intersection_gene



# get markers only for cluster 6
cluster6.markers <- FindMarkers(seurat.wnn.peak,assay = 'peaks_BC7_h3k27ac', ident.1 = 6, logfc.threshold = 0.25, test.use = "wilcox")
# head(cluster2.markers, n = 5)


# plot markers
features_spatial <- c(intersection_gene[1:10], 'Hand2')
features_spatial
DefaultAssay(object = seurat.wnn.peak) <- "gscore_BC7_h3k27ac"

# for multiple feature
for (feature in features_spatial){
  p <- SpatialPlot_new(seurat.wnn.peak, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  print(p)
}

plot_features <- function(feature){
  p <- SpatialPlot_new(seurat.wnn.peak, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_spatial, plot_features)
ggList[[1]]

lapply(seq(length(ggList)),
       function(x)ggsave(filename=paste0('./markers_plot/gscore_BC7_h3k27ac_', features_spatial[x],".png"), plot=ggList[[x]]))

