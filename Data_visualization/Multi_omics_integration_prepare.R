####### This code is for preparing the integrated Seurat object for multi-omics analysis

# basic settings
samplepath_list <- c('/mnt/d/Users/Liran/32nanobody/Spnanob32_H3K27me3/Save-inTissue-Spnanob32_H3K27me3',
                     '/mnt/d/Users/Liran/32nanobody/Spnanob32_H3K4me3/Save-inTissue-Spnanob32_H3K4me3',
                     '/mnt/d/Users/Liran/32nanobody/Spnanob32_ATAC/Save-inTissue-Spnanob32_ATAC')
modalitylist <- c('Spnanob32_H3K27me3','Spnanob32_H3K4me3', 'Spnanob32_ATAC')
colnameslist<-c('Spnanob32_H3K27me3#','Spnanob32_H3K4me3#', 'Spnanob32_ATAC#')
outputname <-'Spananob32_deep_combine'
seurat.rna <- readRDS(file = "/mnt/d/Users/Liran/32nanobody/Spnanob32_RNA/seurat.rna.rds")
seurat.wnn.peak <- seurat.rna

# add the peak matrix and gene score matrix
for (i in 1:3){
  samplepath <- samplepath_list[i]
  modality <- modalitylist[i]
  colname <- colnameslist[i]
  projCUTA <- loadArchRProject(path = samplepath, force = FALSE, showLogo = TRUE)
  
  # call peaks
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
  
  # create a new assay to store ADT information
  adt_assay <- CreateAssayObject(counts = SummarizedExperiment_assay)
  # add this assay to the previously created Seurat object
  seurat.wnn.peak[[paste0('peaks_',modality)]] <- adt_assay
}

# Do intersection with RNA and add gene score 
gene_rna <- rownames(seurat.wnn.peak.rna@assays$Spatial)

for (path in samplepath_list){
  projCUTA <- loadArchRProject(path = path, force = FALSE, showLogo = TRUE)
  projCUTA
  markerGenes <- getFeatures(ArchRProj = projCUTA, useMatrix = "GeneScoreMatrix")
  gene_final <- intersect(gene_rna, markerGenes)
}

for (i in 1:3){
  print(i)
  path <- samplepath_list[i]
  sampleNames <- modalitylist[i]
  gene_score2 <- readRDS(paste0('../Spnano32_combine/', sampleNames,'_all_gene_by_intersection_gene_score.rds'))
  adt_assay <- CreateAssayObject(counts = gene_score2)
  seurat.wnn.peak[[paste0('gscore_intersect_',modality)]] <- adt_assay
  print(seurat.wnn.peak[[paste0('gscore_intersect_',modality)]][gene_to_look, ])
}

saveRDS(seurat.wnn.peak, file = "./seurat.wnn.peak.rna.rds")


