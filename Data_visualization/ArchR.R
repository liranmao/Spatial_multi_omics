####### This code is for basic analysis on ATAC/histone modifications

library(ArchR)
library(Seurat)
library(grid)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(dplyr)

source('./getGeneScore_ArchR.R')
source('./SpatialPlot_new.R')

########## archr project creation
threads = 8
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- './Spnanob38.fragments.sort.bed.gz' # put the fragments file here 
sampleNames <- 'Spananob38_deep_H3K27ac' 

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0, #Dont set this too high because you can always increase later
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme1

saveArchRProject(ArchRProj = projHeme1, outputDirectory = paste0("Save-", sampleNames), load = FALSE)
projHeme1 <- loadArchRProject(path = paste0("Save-", sampleNames), force = FALSE, showLogo = TRUE)


## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = projHeme1))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)

## export meta data
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
summary(meta.data.spatial$nFrags)
summary(meta.data.spatial$Proportion_of_TSS_fragments)
summary(meta.data.spatial$Proportion_of_mito_fragments)
summary(meta.data.spatial$FRiP)

write.table(meta.data.spatial, paste0(sampleNames, '_meta_data.txt'), row.names = TRUE, col.names = TRUE)

## filter off-tissue tixels using image data
projCUTA <- projHeme1[meta.data.spatial$cellID_archr, ]
projCUTA

p <- plotFragmentSizes(ArchRProj = projCUTA)
plotPDF(p, name = paste0(sampleNames, '_FragmentSizes.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)


## dimentsion reduction, clustering, and add umap emedding
projCUTA <- addIterativeLSI(
  ArchRProj = projCUTA,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

projCUTA <- addClusters(
  input = projCUTA,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, # default:0.5
  force = TRUE
)

projCUTA <- addUMAP(
  ArchRProj = projCUTA, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

projCUTA <- addImputeWeights(projCUTA)


saveArchRProject(ArchRProj = projCUTA, outputDirectory = paste0("Save-inTissue-", sampleNames), load = FALSE)
projCUTA <- loadArchRProject(path = paste0("Save-inTissue-", sampleNames), force = FALSE, showLogo = TRUE)
projCUTA

## QC plot
df <- getCellColData(projCUTA, select = c("log10(nFrags)", "TSSEnrichment"))
df

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(100), 4.5),
  ylim = c(1, 8),
  baseSize = 12
)
p
dev.off()

# this will plot under the Save-inTissue- folder, make sure that you load from the path first
plotPDF(p, name = paste0(sampleNames, "TSS-vs-Frags.pdf"), ArchRProj = projCUTA, addDOC = FALSE)


## get markers
markersGS <- getMarkerFeatures(
  ArchRProj = projCUTA, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR < 0.05 & Log2FC >= 0.5")
markerList_pos$C1

for (i in seq_len(length(markerList_pos))) {
  write.table(markerList_pos[[i]], file=paste0('./markers_list/', sampleNames, '_C', i, '_markers.txt'),
              quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
}

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC <= -1",
  transpose = TRUE,
  # invert = FALSE
)

heatmapGS

heatmapGS_complex<-ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS_complex, name = paste0(sampleNames, "matker_motif_heatmap.pdf"), ArchRProj = projCUTA, addDOC = FALSE)


markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)


## get gene score matrix
meta.data <- as.data.frame(getCellColData(ArchRProj = projCUTA))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

projCUTA <- addImputeWeights(projCUTA)
gene_score <- getGeneScore_ArchR(ArchRProj = projCUTA, name = markerGenes, imputeWeights = getImputeWeights(projCUTA))

saveRDS(gene_score, paste0(sampleNames,'_gene_score.rds'))

dim(gene_score)


## For spatial plot visualization
## create seurat object
object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

p2 <- VlnPlot(spatial.obj, features = "nFrags", pt.size = 0.1, log = TRUE) + NoLegend()
median(meta.data$nFrags)
plotPDF(p2, name = paste0(sampleNames, '_nfrags.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)


p <- SpatialPlot_new(spatial.obj, features = "nFrags",  pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
plotPDF(p, name = paste0(sampleNames, '_nFrags_spatial.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)
dev.off()


## cluster and plot
n_clusters <- length(unique(projCUTA$Clusters))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0('C', seq_len(n_clusters))
cols
p1 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4.5, cols = cols, image.alpha = 0, stroke = 0)
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)

plotPDF(p1, name = paste0(sampleNames, '_clusters_spatial.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)



p2 <- plotEmbedding(ArchRProj = projCUTA, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 0.5, pal = cols)

plotPDF(p2, name = paste0(sampleNames, '_clusters_umap.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)


## plot marker features
# for one feature
features_spatial <- c('Wwp2', 'Rhcg', 'Rbfox3', 'Slc4a1',
                      # 'Cdh15', 'Tmem145',
                      'Epb42')
feature <- features_spatial[1]

p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

png(filename = paste0('./markers_plot/', feature, '_upReg_spatial.png'), width = 1200, height = 1200, res = 300)
p
dev.off()

# for multiple feature
plot_features <- function(feature){
  p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}
ggList <- lapply(features_spatial, plot_features)
ggList[[1]]

lapply(seq(length(ggList)),
       function(x)ggsave(filename=paste0('./markers_plot/', features_spatial[x],".png"), plot=ggList[[x]]))



# plot track
p_track <- plotBrowserTrack(
  ArchRProj = projCUTA,
  groupBy = "Clusters",
  geneSymbol = features_spatial,
  upstream = 50000, #50000,
  downstream = 50000, #50000,
  normMethod = 'ReadsInTSS',
  minCells = 0,
  borderWidth = 0
)
plotPDF(plotList = p_track,
        name = paste0(sampleNames, "_Plot-Tracks-Marker-Genes-10kb.pdf"),
        ArchRProj = projCUTA,
        addDOC = FALSE, width = 5, height = 5)






## Integration with scRNA-seq
MOCA_dir <- "/mnt/e/Users/Pengfei/Processed_data/2023/Nanobody_project/Spnb17/Spnanob17_h3k27ac/fragments/integrate_scRNA/"

meta.data.RNA <- read.csv(file = paste0(MOCA_dir, 'cell_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- read.csv(file = paste0(MOCA_dir, 'gene_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- gene.ANN.RNA[, 'gene_short_name', drop = FALSE]

cds <- readRDS(paste0(MOCA_dir, 'gene_count_cleaned_sampled_100k.RDS'))

MOCA <- CreateSeuratObject(counts = cds, project = 'MOCA')
meta.data.RNA <- meta.data.RNA[colnames(MOCA), ]
meta.data.RNA <- meta.data.RNA[, c('Main_cell_type', 'development_stage')]

MOCA <- AddMetaData(object = MOCA, metadata = meta.data.RNA)
MOCA_E11 <- subset(MOCA, development_stage == 13.5) # change the day
MOCA_E11.raw.data <- as.matrix(GetAssayData(MOCA_E11, slot = 'counts'))
MOCA_E11.raw.data <- as.data.frame(MOCA_E11.raw.data)
MOCA_E11.raw.data <- merge(gene.ANN.RNA, MOCA_E11.raw.data, by=0, all=TRUE)
which(is.na(MOCA_E11.raw.data$gene_short_name))

# Remove duplicated names, only keep the first one
tt <- table(MOCA_E11.raw.data$gene_short_name)
name_rep <- names(which(tt > 1))
row_del_fun <- function(x){
  rows <- which(MOCA_E11.raw.data$gene_short_name == x)
  return(rows[2:length(rows)] )
}
row_del <- unlist(lapply(name_rep, row_del_fun))
MOCA_E11.raw.data <- MOCA_E11.raw.data[-row_del, ]
#
row.names(MOCA_E11.raw.data) <- MOCA_E11.raw.data$gene_short_name
MOCA_E11.raw.data <- MOCA_E11.raw.data[, -c(1:2), drop=FALSE]
MOCA_E11 <- CreateSeuratObject(counts = MOCA_E11.raw.data, project = 'MOCA_E11', meta.data = MOCA_E11@meta.data)

projCUTA <- addGeneIntegrationMatrix(
  ArchRProj = projCUTA, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = MOCA_E11,
  addToArrow = TRUE,
  groupRNA = "Main_cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)

# export meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = projCUTA))
meta.data
#write.table(meta.data, paste0(sampleNames, '_meta_data_scRNA_integration.txt'), row.names = TRUE, col.names = TRUE)
#

# pal <- paletteDiscrete(values = MOCA_E11$Main_cell_type)
# pal
# p1 <- plotEmbedding(
#   projCUTA, 
#   colorBy = "cellColData", 
#   name = "predictedGroup_Un", 
#   pal = pal
# )
# p1
# png(filename = 'MOCA_label_UMAP.png', width = 3600, height = 2400, res = 300)
# p1
# dev.off()

meta.data.integration <- as.data.frame(getCellColData(ArchRProj = projCUTA))[, c('predictedCell_Un', 'predictedGroup_Un', 'predictedScore_Un')]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

p <- SpatialDimPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'predictedGroup_Un', pt.size.factor = 4)
p
# png(filename = 'MOCA_label_spatial.png', width = 3600, height = 3000, res = 300)
# p
# dev.off()
library(ggplot2)
library(patchwork)
source('./scripts/SpatialDimPlot_new.R')
Idents(spatial.obj) <- 'predictedGroup_Un'

table(spatial.obj$predictedGroup_Un)
ids.highlight <- names(table(spatial.obj$predictedGroup_Un))[15]
ids.highlight
# p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = ids.highlight), 
#                         facet.highlight = TRUE, pt.size.factor = 3, alpha = c(1,0), stroke = 0)
# p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
# p
#
p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = ids.highlight),
                        facet.highlight = TRUE, pt.size.factor = 3, alpha = c(1,0.05), stroke = 0, image.alpha = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

png(filename = paste0('MOCA_label_highlight_', ids.highlight, '.png'), width = 3600, height = 3600, res = 300)
p
dev.off()

plotPDF(p, name = paste0('MOCA_label_highlight_', ids.highlight, '.pdf'), ArchRProj = projCUTA, addDOC = FALSE)



# plot list
features_spatial <- names(table(spatial.obj$predictedGroup_Un))

# plot_features <- function(feature){
#   p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = feature), 
#                           facet.highlight = TRUE, pt.size.factor = 3, alpha = c(1,0), stroke = 0)
#   p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
#   p
# }

plot_features <- function(feature){
  p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = feature), 
                          facet.highlight = TRUE, pt.size.factor = 3, alpha = c(1,0.05), stroke = 0, image.alpha = 0)
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_spatial, plot_features)
ggList[[1]]

# lapply(seq(length(ggList)), 
#        function(x)ggsave(filename=paste0('./Integration_with_scRNA-seq/', features_spatial[x],".png"), plot=ggList[[x]]))

plotPDF(ggList, name = 'MOCA_label_highlight.pdf', ArchRProj = projCUTA, addDOC = FALSE)
#


## plotTrajectory
source('./scripts/SpatialPlot_traj.R')
table(spatial.obj$predictedGroup_Un)

trajectory <- c("Radial glia", "Postmitotic premature neurons", "Excitatory neurons") #c("Radial glia", "Excitatory neurons")
trajectory

projCUTA <- addTrajectory(
  ArchRProj = projCUTA, 
  name = "Neuron_U", 
  groupBy = "predictedGroup_Un",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)

p <- plotTrajectory(projCUTA, trajectory = "Neuron_U", colorBy = "cellColData", name = "Neuron_U")
p[[1]]

trajMM  <- getTrajectory(ArchRProj = projCUTA, name = "Neuron_U", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p1
plotPDF(p1, name = "Traj_Neuron-Heatmaps_05272021.pdf", ArchRProj = projCUTA, addDOC = FALSE, width = 6, height = 8)


meta.data.integration <- as.data.frame(getCellColData(ArchRProj = projCUTA))[, c('Neuron_U'), drop=FALSE]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names
all(row.names(meta.data.integration) == colnames(spatial.obj))

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)


p <- SpatialPlot_traj(spatial.obj, features = "Neuron_U",  pt.size.factor = 4, image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

png(filename = paste0('Traj_Neuron_05302021.png'), width = 3600, height = 3600, res = 300)
p
dev.off()
library(hexbin)
markers_neuron <- c('Sox2', 'Pax6', 'Ntng1', 'Car10')
markers_neuron[1]
p1 <- plotTrajectory(projCUTA, trajectory = "Neuron_U", colorBy = "GeneScoreMatrix", name = "Sox2", continuousSet = "horizonExtra")
p1[[2]]

# plot list
plot_features <- function(feature){
  p1 <- plotTrajectory(projCUTA, trajectory = "Neuron_U", colorBy = "GeneScoreMatrix", name = feature, continuousSet = "horizonExtra")
  p1[[2]]
}

ggList <- lapply(markers_neuron, plot_features)
ggList[[1]]

plotPDF(ggList, name = 'Traj_genes_05272021.pdf', ArchRProj = projCUTA, addDOC = FALSE)
#


# Co-embedding with scRNA-seq
#source('getGeneIntegrationMatrix__ArchR.R')

MOCA_E11 <- NormalizeData(MOCA_E11)
MOCA_E11 <- FindVariableFeatures(object = MOCA_E11)
MOCA_E11 <- ScaleData(MOCA_E11)

genes.use <- VariableFeatures(MOCA_E11)

getAvailableMatrices(projCUTA)
genes_CUTA <- getFeatures(
  ArchRProj = projCUTA,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)

genes.use <- genes.use[which(genes.use %in% genes_CUTA)]

projCUTA <- addImputeWeights(projCUTA)
imputation <- getGeneScore_ArchR(ArchRProj = projCUTA, name = genes.use, imputeWeights = getImputeWeights(projCUTA)) #name = gsub("\\-", "_", genes.use)
imputation <- log(imputation + 1) #use natural log

object_CUTA <- Seurat::CreateSeuratObject(counts = imputation[head(seq_len(nrow(imputation)), 5), , drop = FALSE], meta.data = meta.data.integration, project = 'CUTA')
object_CUTA[["GeneScore"]] <- Seurat::CreateAssayObject(counts = imputation)

#object_CUTA <- CreateSeuratObject(counts = imputation, assay = 'GeneScore', meta.data = meta.data.integration, project = 'CUTA')
DefaultAssay(object_CUTA) <- "GeneScore"
#object_CUTA <- NormalizeData(object_CUTA)
object_CUTA <- Seurat::ScaleData(object_CUTA)

transfer.anchors <- FindTransferAnchors(reference = MOCA_E11, query = object_CUTA, features = VariableFeatures(object = MOCA_E11), 
                                        reference.assay = "RNA", query.assay = "GeneScore", reduction = "cca")

rD <- getReducedDims(ArchRProj = projCUTA, reducedDims = "IterativeLSI", corCutOff = 0.75, dimsToUse = 1:30)
new_row_names <- row.names(rD)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(rD) <- new_row_names

weight.reduction <- CreateDimReducObject(
  embeddings = rD, 
  key = "LSI_", 
  assay = DefaultAssay(object_CUTA)
)

refdata <- GetAssayData(MOCA_E11, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = weight.reduction, dims = seq_len(ncol(rD)))

object_CUTA[["RNA"]] <- imputation

coembed <- merge(x = MOCA_E11, y = object_CUTA)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)

coembed <- FindNeighbors(coembed, dims = 1:30)
coembed <- FindClusters(coembed, resolution = 0.1)

coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$Main_cell_type), coembed$Main_cell_type, coembed$predictedGroup_Un)
coembed$CUTA_cells <- ifelse(!is.na(coembed$Main_cell_type), coembed$Main_cell_type, 'zCUTA')

coembed$MOCA_cells <- ifelse(is.na(coembed$Main_cell_type), coembed$predictedGroup_Un, 'zMOCA')

n_clusters <- length(unique(coembed$CUTA_cells))
cols <- scales::hue_pal()(n_clusters)
names(cols) <- unique(coembed$CUTA_cells)
cols
cols['zCUTA']
cols['zCUTA'] <- '#000000'
cols

# p <- DimPlot(coembed, group.by = "orig.ident", cols = c('CUTA'='black'))
# p <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE, split.by = 'orig.ident')
# p <- DimPlot(coembed, label = TRUE, repel = TRUE, split.by = 'orig.ident')
p <- DimPlot(coembed, group.by = "CUTA_cells", cols = cols, pt.size = 1.2)
p

png(filename = 'MOCA_coembed_UMAP.png', width = 4500, height = 3000, res = 300)
p
dev.off()

#n_clusters <- length(unique(coembed$MOCA_cells))
#cols <- scales::hue_pal()(n_clusters)
names(cols)[38] <- 'zMOCA'
cols
cols['zMOCA']
cols['zMOCA'] <- '#BFBFBF'
cols

p <- DimPlot(coembed, group.by = "MOCA_cells", cols = cols, pt.size = 1.2)
p

png(filename = 'MOCA_coembed_UMAP_MOCA_cells.png', width = 4500, height = 3000, res = 300)
p
dev.off()

