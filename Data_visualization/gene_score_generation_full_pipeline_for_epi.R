rm(list=ls())
library(ArchR)
library(Seurat)
library(grid)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(dplyr)

########## Setup paths
base_dir <-   # <-- CHANGE THIS: where your sample folders are
output_dir <-          # <-- CHANGE THIS: where to save results
script_dir <- '/home/liran/processed_data/Users/Liran/Processed_data/2025/02_nano_review/17nanobody_deep/Spnanob17_h3k27me3_deep/scripts'        # <-- CHANGE THIS: where your helper scripts are

source(file.path(script_dir, 'getGeneScore_ArchR.R'))
source(file.path(script_dir, 'SpatialPlot_new.R'))

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

########## ArchR settings
threads <- 8
addArchRThreads(threads = threads)
addArchRGenome("mm10")

########## Define all samples
samples <- c("P5S2", "P21S2")
samples <- c("P7S2")

########## Function to process one sample
process_sample <- function(sample_name, base_dir, output_dir) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Processing sample:", sample_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  # Define paths
  sample_dir <- file.path(base_dir, sample_name, "ATAC")
  spatial_dir <- file.path(sample_dir, "spatial")
  fragments_file <- file.path(sample_dir, "atac_fragments.tsv.gz")
  
  # Check if files exist
  if (!file.exists(fragments_file)) {
    cat("WARNING: Fragment file not found for", sample_name, "- Skipping\n")
    return(NULL)
  }
  
  if (!dir.exists(spatial_dir)) {
    cat("WARNING: Spatial directory not found for", sample_name, "- Skipping\n")
    return(NULL)
  }
  
  # Set working directory to sample ATAC folder
  setwd(sample_dir)
  
  ########## Create Arrow Files
  ArrowFiles <- createArrowFiles(
    inputFiles = fragments_file,
    sampleNames = sample_name,
    filterTSS = 0,
    filterFrags = 0,
    minFrags = 0,
    maxFrags = 1e+07,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    offsetPlus = 0,
    offsetMinus = 0,
    TileMatParams = list(tileSize = 5000)
  )
  
  ########## Create ArchR Project
  projHeme1 <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = file.path(sample_dir, paste0(sample_name, "_ArchRProject")),
    copyArrows = TRUE
  )
  
  ########## Prepare meta data
  meta.data <- as.data.frame(getCellColData(ArchRProj = projHeme1))
  meta.data['cellID_archr'] <- row.names(meta.data)
  new_row_names <- row.names(meta.data)
  new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
  new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
  row.names(meta.data) <- new_row_names
  
  ########## Read spatial image and filter
  assay <- "Spatial"
  filter.matrix <- TRUE
  slice <- "slice1"
  
  image <- Read10X_Image(image.dir = spatial_dir, filter.matrix = filter.matrix)
  
  # Filter to on-tissue spots
  on_tissue_barcodes <- intersect(row.names(meta.data), Cells(image))
  meta.data.spatial <- meta.data[on_tissue_barcodes, ]
  
  cat("Total barcodes:", nrow(meta.data), "\n")
  cat("On-tissue barcodes:", length(on_tissue_barcodes), "\n")
  
  # Filter ArchR project
  projCUTA <- projHeme1[meta.data.spatial$cellID_archr, ]
  
  ########## Dimension reduction and clustering
  projCUTA <- addIterativeLSI(
    ArchRProj = projCUTA,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list(
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
    resolution = 1,
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
  
  ########## Get gene scores
  meta.data <- as.data.frame(getCellColData(ArchRProj = projCUTA))
  meta.data['cellID_archr'] <- row.names(meta.data)
  new_row_names <- row.names(meta.data)
  new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
  new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
  row.names(meta.data) <- new_row_names
  
  gene_all <- getFeatures(ArchRProj = projCUTA, useMatrix = "GeneScoreMatrix")
  gene_score <- getGeneScore_ArchR(
    ArchRProj = projCUTA, 
    name = gene_all, 
    imputeWeights = getImputeWeights(projCUTA)
  )
  
  # Save gene score matrix
  saveRDS(gene_score, file.path(output_dir, paste0(sample_name, '_gene_score_all.rds')))
  
  ########## Create Seurat object
  object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)
  
  image <- Read10X_Image(image.dir = spatial_dir, filter.matrix = filter.matrix)
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  
  spatial.obj <- object
  
  ########## Generate plots
  # Cluster colors
  n_clusters <- length(unique(projCUTA$Clusters))
  cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
  names(cols) <- paste0('C', seq_len(n_clusters))
  
  # Spatial plot
  p1 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, 
                    group.by = 'Clusters', pt.size.factor = 1, 
                    cols = cols, image.alpha = 1, stroke = 0)
  p1$layers[[1]]$aes_params$shape <- 22
  
  png(file.path(output_dir, paste0(sample_name, '_spatial_plot.png')), 
      width = 800, height = 600)
  print(p1)
  dev.off()
  
  # UMAP plot
  p2 <- plotEmbedding(ArchRProj = projCUTA, colorBy = "cellColData", 
                      name = "Clusters", embedding = "UMAP", size = 0.5)
  
  png(file.path(output_dir, paste0(sample_name, '_dim_plot.png')), 
      width = 800, height = 600)
  print(p2)
  dev.off()
  
  # Save Seurat object
  saveRDS(spatial.obj, file.path(output_dir, paste0(sample_name, '_seurat_obj_all.rds')))
  
  # Save ArchR project
  saveArchRProject(ArchRProj = projCUTA, 
                   outputDirectory = file.path(output_dir, paste0(sample_name, "_ArchRProject")),
                   load = FALSE)
  
  cat("\nCompleted processing:", sample_name, "\n")
  
  return(spatial.obj)
}

########## Run for all samples
results <- list()

for (sample in samples) {
  tryCatch({
    results[[sample]] <- process_sample(sample, base_dir, output_dir)
  }, error = function(e) {
    cat("\nERROR processing", sample, ":", conditionMessage(e), "\n")
    results[[sample]] <- NULL
  })
}

########## Summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PROCESSING COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

successful <- names(results)[!sapply(results, is.null)]
failed <- names(results)[sapply(results, is.null)]

cat("\nSuccessfully processed:", length(successful), "samples\n")
if (length(successful) > 0) cat("  ", paste(successful, collapse = ", "), "\n")

if (length(failed) > 0) {
  cat("\nFailed samples:", length(failed), "\n")
  cat("  ", paste(failed, collapse = ", "), "\n")
}

cat("\nOutput saved to:", output_dir, "\n")

