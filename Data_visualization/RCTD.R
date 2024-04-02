####### This code is for using RTCD for cell type deconvolution

library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# extract information to pass to the RCTD Reference function
ref <- MOCA_E11_filtered
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$Main_cell_type)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster)#, nUMI)

# set up query with the RCTD function SpatialRNA
slide.seq <- spatial.obj
counts <- slide.seq[["Spatial"]]$counts
coords <- GetTissueCoordinates(slide.seq)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))


RCTD <- create.RCTD(query, reference, max_cores = 4)#, MAX_MULTI_TYPES = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

# Create the output directory in your working directory
resultsdir <- "RCTD_output_plots/" # An output directory for plots, can be any name
dir.create(resultsdir)

# Create variables from the myRCTD object to plot results
barcodes <- colnames(RCTD@spatialRNA@counts) # list of spatial barcodes
weights <- RCTD@results$weights # Weights for each cell type per barcode

# Normalize per spot weights so cell type probabilities sum to 1 for each spot
norm_weights <- normalize_weights(weights)
cell_type_names<-colnames(norm_weights) # List of cell types

# Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
# Save each plot as a jpg file
for(i in 1:length(cluster)){
  plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=1)
  ggsave(paste(resultsdir, cell_type_names[i],'_weights.jpg', sep=''), height=5, width=5, units='in', dpi=300)
}

deconv_est <- as.matrix(norm_weights)
