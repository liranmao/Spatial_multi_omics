####### This code is for using FigR for gene regulation analysis

library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ComplexHeatmap)
library(networkD3)


## Preparation
# Extract counts matrix and features from ATAC
counts_matrix <- seurat.wnn@assays$peaks_H3K27ac$counts
features <- seurat.wnn@assays$peaks_H3K27ac@ranges

# Create SummarizedExperiment object
ATAC.se <- SummarizedExperiment(assays=list(counts=counts_matrix), rowRanges=features)
dim(ATAC.se) # Peaks x Cells

# RNA
DefaultAssay(seurat.wnn) <- 'RNA'
seurat.wnn <- NormalizeData(seurat.wnn)
RNAmat <- seurat.wnn@assays$RNA@data
dim(RNAmat) # Genes x Cells

RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]
dim(RNAmat) # Genes x Cells

## FigR
# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                                 RNAmat = RNAmat,
                                 genome = "mm10", # One of hg19, mm10 or hg38 
                                 nCores = 30,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = 100)
saveRDS(cisCorr,file = 'cisCorr_03182024_FigR.rds')
head(cisCorr)


cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff = 5, # No. sig peaks needed to be called a DORC, cutoff=2
                       labelTop = 25,
                       returnGeneList = TRUE, # Set this to FALSE for just the plot
                       labelSize = 2,
                       force=2)

# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

numDorcs_cutoff5 <- numDorcs[which(numDorcs$n>=5),]

cellkNN <- seurat.wnn@neighbors$weighted.nn@nn.idx
dim(cellkNN)
rownames(cellkNN) <- seurat.wnn@neighbors$weighted.nn@cell.names

dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 30)
dim(dorcMat)

# Smooth dorc scores using cell KNNs (k=30)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN,mat = dorcMat,nCores = 10)
saveRDS(dorcMat.s, file = 'dorcMat.s.rds')

# Smooth RNA using cell KNNs
RNAmat.s <- smoothScoresNN(NNmat = cellkNN,mat = RNAmat,nCores = 4)
saveRDS(RNAmat.s, file = 'RNAmat.s.rds')

# Visualize on pre-computed UMAP
umap.d <- as.data.frame(Embeddings(seurat.wnn, reduction = "wnn.umap"))

# DORC score for Col2a1
dorcg <- plotMarker2D(umap.d,dorcMat.s,markers = c("Zic4"),maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("Zic4 DORC")
dorcg


figR.d <- runFigRGRN(ATAC.se = ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "mm10",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s, 
                     nCores = 20)
saveRDS(figR.d,file = 'figR.d_03182024_FigR_cutoff5.rds')


# figR.d <- readRDS('./figR.d_03142024_FigR.rds')
figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

rankDrivers(figR.d,rankBy = "meanScore",interactive = FALSE)
rankDrivers(figR.d,score.cut = 1.5,rankBy = "nTargets",interactive = TRUE)
plotDrivers(figR.d,score.cut = 1,marker = "Ap2a1")


plotfigRHeatmap(figR.d = figR.d,
                score.cut = 1,
                TFs = c("Lef1","Dlx3","Grhl1","Gata6","Klf3","Barx2","Pou2f3"),
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
)





p1 <- plotfigRNetwork(figR.d,
                      score.cut = 1.5,
                      TFs = c('Neurod2','Mga'),
                      weight.edges = TRUE)
p1

p2 <- plotfigRNetwork(figR.d,
                      score.cut = 1.5,
                      TFs = c('Hoxc5'),
                      weight.edges = TRUE)
p2


## Spatial visualization
# create a new assay to store DORC information
dorc_assay <- CreateAssay5Object(counts = dorcMat.s)

# add this assay to the previously created Seurat object
spatial_RNA[["DORC"]] <- dorc_assay

DefaultAssay(spatial_RNA) <- 'DORC'

dorc_features <- c('Sox10') 
#dorc_features <- c('Lhx5', 'Nes', 'Neurod2', 'Nrxn2')
dorc_features <- numDorcs_cutoff5$Gene


p <- SpatialFeaturePlot(object = spatial_RNA, slot = 'counts', features = dorc_features, alpha = c(0.5, 1), ncol = 1, 
                        pt.size.factor = 4.5, image.alpha = 0, min.cutoff = "q10", max.cutoff = "q90", stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

png(filename = paste0('./FigR_plots/', dorc_features, '_spatial.png'), width = 2400, height = 2400, res = 300)
p
dev.off()


# plot list
features_spatial <- dorc_features

plot_features <- function(feature){
  p <- SpatialFeaturePlot(object = spatial_RNA, slot = 'counts', features = feature, alpha = c(0.5, 1),
                          pt.size.factor = 4.5, image.alpha = 0, min.cutoff = "q10", max.cutoff = "q90", stroke = 0)
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_spatial, plot_features)
ggList[[1]]

lapply(seq(length(ggList)),
       function(x)ggsave(filename=paste0('./figR/dorc_plot/', features_spatial[x],".png"), plot=ggList[[x]]))


# visualization for corresponding RNA expression
library(Nebulosa)
DefaultAssay(spatial_RNA) <- 'Spatial_RNA'
SpatialFeaturePlot(object = spatial.obj, features = dorc_features, ncol = 3, 
                   pt.size.factor = 4.5, image.alpha = 0, stroke = 0)

plot_density(spatial.obj, 'Neurod2')


DefaultAssay(seurat.wnn) <- 'RNA'

pd <- plot_density(seurat.wnn, "Dcx", reduction = 'wnn.umap')
pd

pd_data <- pd$data
pd_data <- pd_data[, c(3), drop=FALSE]

all(Cells(spatial.obj) == row.names(pd_data))
pd_data <- pd_data[Cells(spatial.obj), , drop=FALSE]
all(Cells(spatial.obj) == row.names(pd_data))

spatial.obj <- AddMetaData(object = spatial.obj, metadata = pd_data)


p <- SpatialFeaturePlot(spatial.obj, features = "feature",  pt.size.factor = 4.5, image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

