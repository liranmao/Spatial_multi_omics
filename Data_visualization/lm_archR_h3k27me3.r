########## tutorial
# https://docs.atlasxomics.com/projects/AtlasXpipeline/en/latest/SpatialEpigenome.html

rm(list=ls())

########## prepare and load packages
# Installation of ArchR
# install.packages("devtools")
# 
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
# ArchR::installExtraPackages()
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# devtools::install_github("GreenleafLab/chromVARmotifs")

library(ArchR)
library(Seurat)
library(grid)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(dplyr)


########## creat noe folders to store the data
getwd()
setwd('/mnt/d/Users/Liran/17nanobody/Spnanob17_H3K27me3/outs')
# Specify the path to the new folder/directory, markers_list, markers_go, markers_list
for (i in c('markers_list', 'markers_go','maker_peaks','markers_plot','motifs_list','motifs_plot')){
  folder_path <- paste0("./",i)
  # Create the new folder
  dir.create(folder_path)
  
}

source('./scripts/getGeneScore_ArchR.R')
source('./scripts/SpatialPlot_new.R')

########## archr project creation
threads = 8
addArchRThreads(threads = threads)

addArchRGenome("mm10")

# bc5 is atac
inputFiles <- './fragments.tsv.gz'
sampleNames <- 'Spananob17_BC6'

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
# projHeme1 <- loadArchRProject(path = paste0("Save-", sampleNames), force = FALSE, showLogo = TRUE)


############ spatial object creation
## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = projHeme1))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

qc_10x <- read.csv('./singlecell.csv', header = TRUE)
qc_10x['Proportion_of_TSS_fragments'] <- qc_10x$TSS_fragments / qc_10x$passed_filters
qc_10x['Proportion_of_mito_fragments'] <- qc_10x$mitochondrial / qc_10x$total
qc_10x['FRiP'] <- qc_10x$peak_region_fragments / qc_10x$passed_filters

summary(qc_10x$Proportion_of_TSS_fragments)
summary(qc_10x$Proportion_of_mito_fragments)
summary(qc_10x$FRiP)

qc_10x <- qc_10x[-1,]
row.names(qc_10x) <- qc_10x$barcode
new_row_names <- row.names(qc_10x)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(qc_10x) <- new_row_names
#qc_10x <- qc_10x[, c('Proportion_of_TSS_fragments', 'Proportion_of_mito_fragments'), drop =FALSE]
qc_10x <- qc_10x[, c('Proportion_of_TSS_fragments', 'Proportion_of_mito_fragments', 'FRiP'), drop =FALSE]

meta.data <- merge(meta.data, qc_10x, by=0, all=TRUE)
row.names(meta.data) <- meta.data$Row.names
meta.data <- meta.data[-1]

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

# png(filename = 'lm_FragmentSizes.png', width = 1200, height = 1200, res = 300)
# p
dev.off()


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

# heatmapGS <- markerHeatmap(
#   seMarker = markersGS,
#   cutOff = "FDR <= 0.05 & Log2FC <= -1",
#   transpose = TRUE
# )

# heatmapGS <- markerHeatmap(
#   seMarker = markersGS,
#   cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
#   transpose = TRUE
# )

heatmapGS_complex<-ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS_complex, name = paste0(sampleNames, "matker_motif_heatmap.pdf"), ArchRProj = projCUTA, addDOC = FALSE)


# ?????
# markers_cluster <- as.data.frame(markerList_neg$C4)
# markers_cluster


markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)
#markerGenes <- getFeatures(ArchRProj = projCUTA, useMatrix = "GeneScoreMatrix")

# markerGenes <- c('Nfe2', 'Lyl1', 'Tal1', 'Kcnn4', 'Slc4a1', 'Gfi1b', 'Vamp8', 'Kel',
#                  'Cnp', 'Aatk', 'Sox10', 'Mal', 'Gjc2', 'Gal3st1', 'Cldn11', 'Plekhb1',
#                  'Adcy5', 'Pde7b', 'Pde10a', 'Gng7', 'Pde1b', 'Lingo3', 'Rgs9', 'Dlx6',
#                  'Pax6', 'Nr2e1', 'Prdm16', 'Sox1ot', 'Sall1', 'Sox2ot', 'Sox21', 'Sox2',
#                  'Tbx15', 'Adgra2', 'Tbx2', 'Sox17', 'Fli1', 'Hoxc6', 'Tbx3', 'Ets1',
#                  'Lhx1', 'Ptprn2', 'Pou3f2', 'Rnf165', 'Tubb2b', 'Otx2',
#                  'Tbr1', 'Khdrbs3', 'Nptx1', 'Shank2', 'Nfix', 'Actl10', 'Sstr3', 'Slc17a7') # spatial-CUT&Tag K4me3 markers

# markerGenes <- c('Bnc2', 'Sgk2', 'Grik3', 'Casz1', 'Rbms3', 'Ptprm',
#                  'Sez6l', 'Gabbr2', 'Ephb2', 'Olfm1', 'Myo18b', 'Ncald', 'Kcnj6', 'Grk3',
#                  'Rbms3', 'Rarb', 'Rfx4', 'Myo3b', 'Reln', 'Rxrg', 'Six3os1', 'Kif26b',
#                  'Ets1', 'Mgmt', 'Pard3b', 'Suclg2', 'Mecom', 'Dock2', 'Ptpn14', 'Alk',
#                  'Itpr2', 'Zfp536', 'Pacrg', 'Spata31', 'Pde3a', 'Samd4', 'Creb5', 'Tns3') # spatial-CUT&Tag K27me3 markers

# markerGenes <- c('Bcl11b', 'Cartpt', 'Cux1', 'Cux2', 'Fam19a1', 'Foxp2', 'Grik3',
#                  'Pou3f2', 'Pparg', 'Rora', 'Rorb', 'Sdk2', 'Shb', 'Slc17a6',
#                  'Syt6', 'Thsd7b')


## Spatial plots
## Prepare meta data again for in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = projCUTA))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

projCUTA <- addImputeWeights(projCUTA)
gene_score <- getGeneScore_ArchR(ArchRProj = projCUTA, name = markerGenes, imputeWeights = getImputeWeights(projCUTA))
saveRDS(gene_score, paste0(sampleNames,'_gene_score.rds'))


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


# png(filename = paste0(sampleNames,'_nfrags','.png'), width = 1200, height = 1200, res = 300)
# p2
# dev.off()

p <- SpatialPlot_new(spatial.obj, features = "nFrags",  pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
plotPDF(p, name = paste0(sampleNames, '_nFrags_spatial.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)
dev.off()

# plot(0)
# png(filename = paste0(sampleNames, '_nFrags_spatial.png'), width = 2400, height = 2400, res = 300)
# p
# dev.off()

## cluster and plot
n_clusters <- length(unique(projCUTA$Clusters))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0('C', seq_len(n_clusters))
cols
p1 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4.5, cols = cols, image.alpha = 0, stroke = 0)
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)


plotPDF(p1, name = paste0(sampleNames, '_clusters_spatial.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)


# png(filename = paste0(sampleNames, '_clusters_spatial.png'), width = 3600, height = 3000, res = 300)
# p1
# dev.off()

p2 <- plotEmbedding(ArchRProj = projCUTA, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 0.5, pal = cols)

plotPDF(p2, name = paste0(sampleNames, '_clusters_umap.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)


## plot marker features
# for one feature
# features_spatial <- c('Wwp2', 'Rhcg', 'Rbfox3', 'Slc4a1', 
#                       # 'Cdh15', 'Tmem145', 
#                       'Epb42')
# feature <- features_spatial[1]
# 
# p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
#   theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
# p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
# p
# 
# png(filename = paste0('./markers_plot/', feature, '_upReg_spatial.png'), width = 1200, height = 1200, res = 300)
# p
# dev.off()
# 
# # for multiple feature
# plot_features <- function(feature){
#   p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
#     theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
#   p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
#   p
# }
# ggList <- lapply(features_spatial, plot_features)
# ggList[[1]]
# 
# lapply(seq(length(ggList)),
#        function(x)ggsave(filename=paste0('./markers_plot/', features_spatial[x],".png"), plot=ggList[[x]]))



## plot track
# p_track <- plotBrowserTrack(
#   ArchRProj = projCUTA,
#   groupBy = "Clusters",
#   geneSymbol = features_spatial,
#   upstream = 50000, #50000,
#   downstream = 50000, #50000,
#   normMethod = 'ReadsInTSS',
#   minCells = 0,
#   borderWidth = 0
# )
# plotPDF(plotList = p_track,
#         name = paste0(sampleNames, "_Plot-Tracks-Marker-Genes-10kb.pdf"),
#         ArchRProj = projCUTA,
#         addDOC = FALSE, width = 5, height = 5)
# 


## 
