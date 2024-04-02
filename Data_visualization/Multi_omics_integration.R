####### This code is for multi-omics analysis

library(Signac)
library(zellkonverter)
library(Seurat)
library(ggplot2)
library(patchwork)
library(slingshot)
library(ArchR)
library(pheatmap)
library(tidyr)
library(RColorBrewer)


## basic settings
seurat_assay_list <- c('rna_log', 'peaks_Spananob32_H3K27me3', 'peaks_Spnanob32_H3K4me3', 'peaks_Spnanob32_ATAC')
outputname <-'Spananob32_combine'
filter.matrix = TRUE
slice = "slice1"
source('/mnt/d/Users/Liran/32nanobody/Spnanob32_H3K27me3/scripts/getGeneScore_ArchR.R')
source('/mnt/d/Users/Liran/32nanobody/Spnanob32_H3K27me3/scripts/SpatialPlot_new.R')
seurat.wnn.peak <- readRDS("/mnt/d/Users/Liran/32nanobody/Spnano32_combine/seurat.wnn.peak.rna.rds")
seurat.wnn.peak[['rna_raw']] <- seurat.wnn.peak[['Spatial']]


## wnn lsi for peaks and pca for rna
seurat.wnn.peak[['rna_log']] <- seurat.wnn.peak[['Spatial']]
DefaultAssay(seurat.wnn.peak) <- "rna_log"
seurat.wnn.peak <- NormalizeData(seurat.wnn.peak, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.wnn.peak <- FindVariableFeatures(seurat.wnn.peak, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.wnn.peak)
seurat.wnn.peak <- ScaleData(seurat.wnn.peak, features = all.genes)
seurat.wnn.peak <- RunPCA(seurat.wnn.peak, assay = "rna_log", verbose = FALSE, features = VariableFeatures(object = seurat.wnn.peak), reduction.name = 'rna_pca') %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


for(x in c('peaks_Spananob32_H3K27me3', 'peaks_Spnanob32_H3K4me3', 'peaks_Spnanob32_ATAC')){
  DefaultAssay(seurat.wnn.peak) <- x
  seurat.wnn.peak <- RunTFIDF(seurat.wnn.peak) %>% FindTopFeatures() %>% RunSVD(reduction.name = paste0(x,'_lsi'))
}

## examine the dimension components
ElbowPlot(seurat.wnn.peak)
DepthCor(seurat.wnn.peak, assay = 'peaks_Spananob32_H3K27me3', reduction = "peaks_Spananob32_H3K27me3_lsi")
DepthCor(seurat.wnn.peak, assay = 'peaks_Spnanob32_H3K4me3', reduction = "peaks_Spnanob32_H3K4me3_lsi")
DepthCor(seurat.wnn.peak, assay = 'peaks_Spnanob32_ATAC', reduction = "peaks_Spnanob32_ATAC_lsi")


## calculate a WNN graph
seurat.wnn.peak <- FindMultiModalNeighbors(
  seurat.wnn.peak, reduction.list = list('rna_pca', 'peaks_Spananob32_H3K27me3_lsi', 'peaks_Spnanob32_H3K4me3_lsi', 'peaks_Spnanob32_ATAC_lsi'), 
  dims.list = list(1:30,2:30,2:30, 2:30), modality.weight.name = "all.weight",k.nn = 10
)

seurat.wnn.peak <- RunUMAP(seurat.wnn.peak, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

## using different resolution to calculate the clusters
for (reso in c(0.2, 0.8, 1.2, 1.5)){
  seurat.wnn.peak <- FindClusters(seurat.wnn.peak, graph.name = "wsnn", algorithm = 3, resolution = reso, verbose = FALSE)
  n_clusters <- length(unique(seurat.wnn.peak@meta.data[paste0('wsnn_res.', reso)][,1]))
  cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
  names(cols) <- paste0( seq_len(n_clusters)-1)
  
  
  ##### plot 
  # cluster in dimplot
  p3 <- DimPlot(seurat.wnn.peak,label=TRUE,pt.size=0.5, cols = cols,group.by = paste0('wsnn_res.', reso), reduction = 'wnn.umap')
  pdf(paste0("./new_combine_plot/", outputname, "_clusters_umap_wsnn_res.", reso, ".pdf"), width = 5, height = 4)
  print(p3)
  dev.off()
  
  # weight in dimplot
  for(assay in seurat_assay_list){
    assay <- paste0(assay,'.weight')
    p2 <- FeaturePlot(seurat.wnn.peak,assay,min.cutoff = 0.2,max.cutoff = 0.6, reduction = 'wnn.umap') + scale_color_viridis_c()
    pdf(paste0("./new_combine_plot/", outputname, '_', assay,"_umap_wsnn_res.pdf"), width = 5, height = 4)
    print(p2)
    dev.off()
  }
  
  # spatial cluster
  p4 <- SpatialPlot(seurat.wnn.peak, label = FALSE, label.size = 3, group.by = paste0('wsnn_res.', reso), pt.size.factor = 4.5, cols = cols, image.alpha = 0, stroke = 0)
  p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
  p4
  
  pdf(paste0("./new_combine_plot/", outputname, "_clusters_spatial_wsnn_res.", reso, ".pdf"), width = 5, height = 4)
  print(p4)
  dev.off()
  
  # spatial plot weight cluster
  for(assay in seurat_assay_list){
    assay <- paste0(assay,'.weight')
    p3 <- SpatialPlot(seurat.wnn.peak, label = FALSE, label.size = 3, features = assay,pt.size.factor = 4.5,  image.alpha = 0, stroke = 0, min.cutoff = 0.3, max.cutoff = 0.9)+ scale_color_viridis_c()
    p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
    pdf(paste0("./new_combine_plot/", outputname, '_', assay,"_spatial_wsnn_res.pdf"), width = 5, height = 4)
    print(p3)
    dev.off()
  }
}

seurat_assay_list


## lineage analysis
# First, get a subset of seurat object for lineage analysis
seurat.wnn.peak.sub 
group_by_name = 'h3k27me3_cluster'
reduction_fun =  'umap.rna'

# visualize the subset
p3 <- DimPlot(seurat.wnn.peak.sub,label=TRUE,pt.size=0.5,group.by = group_by_name, reduction = reduction_fun)
print(p3)

pdf(paste0("new_combine_plot/", outputname, "_clusters_umap_", group_by_name, '.',reduction_fun, ".pdf"), width = 5, height = 4)
print(p3)
dev.off()

p4 <- SpatialPlot(seurat.wnn.peak.sub, label = FALSE, label.size = 3, group.by = 'wsnn_res.0.2', crop = FALSE, pt.size.factor = 4.5, image.alpha = 0.5, stroke = 0)
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
print(p4)

pdf(paste0("new_combine_plot/", outputname, "_clusters_spatial_", 'wsnn_res.0.2', '.',reduction_fun, ".pdf"), width = 5, height = 5)
print(p4)
dev.off()

# using slingshot
sshot <- slingshot(data = Embeddings(seurat.wnn.peak.sub,reduction = reduction_fun))
pt <- slingPseudotime(sshot)
seurat.wnn.peak.sub@meta.data <- cbind(seurat.wnn.peak.sub@meta.data, pt)

p1 <- FeaturePlot(seurat.wnn.peak.sub,features = 'Lineage1', reduction = reduction_fun) + scale_color_viridis_c()
print(p1)
pdf(paste0("new_combine_plot/", outputname,"_umap_slingshot_", group_by_name, '.',reduction_fun, ".pdf"), width = 5, height = 4)
print(p1)
dev.off()

seurat.wnn.peak.sub$Lineage1[is.na(seurat.wnn.peak.sub$Lineage1)] <- 0
p3 <- SpatialPlot(seurat.wnn.peak.sub, label = FALSE, label.size = 3, features = 'Lineage1',crop = FALSE, pt.size.factor = 4.5,  image.alpha = 0.5, stroke = 0, max.cutoff = 'q98')
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
print(p3)
pdf(paste0("new_combine_plot/", outputname,"_clusters_spatial_slingshot_", group_by_name, '.',reduction_fun, ".pdf"), width = 5, height = 5)
print(p3)
dev.off()


# linear transformation for pt
seurat.wnn.peak.sub$pt <- -seurat.wnn.peak.sub$Lineage1+14
p3 <- SpatialPlot(seurat.wnn.peak.sub, label = FALSE, label.size = 3, features = 'pt',crop = FALSE, pt.size.factor = 4.5,  image.alpha = 0.5, stroke = 0, min.cutoff = 'q2')
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
print(p3)
pdf(paste0("new_combine_plot/", outputname,"_clusters_spatial_slingshot_revers_pt", group_by_name, '.',reduction_fun, ".pdf"), width = 5, height = 5)
print(p3)
dev.off()

DefaultAssay(seurat.wnn.peak.sub) <- 'rna_log'



# dotplot
plot_gene_lineage <- function(genelist, file_name){
  final_list <- genelist
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_H3K27me3@data), final_list)
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_H3K4me3@data), final_list)
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_ATAC@data), final_list)
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$rna_log@data), final_list)  
  
  
  df.plot <- data.frame('cell' = rownames(seurat.wnn.peak.sub@meta.data),
                        'pt' = seurat.wnn.peak.sub$pt,
                        'h3k27me3' = colSums(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_H3K27me3@data[final_list,])/sum(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_H3K27me3@data[final_list,]),
                        'H3K4me3' = colSums(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_H3K4me3@data[final_list,])/sum(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_H3K4me3@data[final_list,]),
                        'ATAC' = colSums(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_ATAC@data[final_list,])/sum(seurat.wnn.peak.sub@assays$gscore2_Spnanob32_ATAC@data[final_list,]))
  
  # 'rna'= colSums(seurat.wnn.peak.sub@assays$rna_log@data[final_list,])/sum(seurat.wnn.peak.sub@assays$rna_log@data[final_list,]))
  df.plot <- df.plot[df.plot$pt>7,]
  df_long <- pivot_longer(
    data = df.plot[-1],
    cols = -pt, # Select all columns except pseudotime to pivot
    names_to = "condition", # This will be the new column for the condition names
    values_to = "score" # This will be the new column for the scores
  )
  
  df_long_2 <- df_long[df_long$pt!=0,]
  
  p <- ggplot(data = df_long_2, aes(x = pt, y = score, color = condition)) +
    geom_point(alpha = 0.5) + # Plots the points with some transparency
    geom_smooth(method = "loess", se = TRUE, aes(fill = condition), alpha = 0.1) + # Adds a smooth line with a confidence interval
    # scale_color_manual(values = c("rna" = "red", "h3k27me3" = "green", "h3k27mac" = "blue")) +
    labs(x = "Pseudo-time", y = "Score", title = "Chromatin opening") +
    theme_minimal() + # Sets a minimal theme
    theme(legend.title = element_blank()) # Removes the legend title
  
  print(p)
  
  pdf(paste0("./new_combine_plot/", outputname,"_clusters_spatial_slingshot_", group_by_name, '.',reduction_fun, "RNA_", file_name, "_dot.pdf"), width = 6.5, height = 7)
  print(p)
  dev.off()
}

length(final_list1)
plot_gene_lineage(final_list1, 'pos_cor_0.1')

plot_one_gene_lineage_no_RNA_nor <- function(gene){
  final_list <- gene
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore_intersect_Spnanob32_H3K27me3@data), final_list)
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore_intersect_Spnanob32_H3K4me3@data), final_list)
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore_intersect_Spnanob32_ATAC@data), final_list)
  
  if (length(final_list) !=0){
    df.plot <- data.frame('cell' = rownames(seurat.wnn.peak.sub@meta.data),
                          'pt' = seurat.wnn.peak.sub$pt,
                          'h3k27me3' = seurat.wnn.peak.sub@assays$gscore_intersect_Spnanob32_H3K27me3@data[final_list,],
                          'h3k4me3' = seurat.wnn.peak.sub@assays$gscore_intersect_Spnanob32_H3K4me3@data[final_list,],
                          'ATAC' = seurat.wnn.peak.sub@assays$gscore_intersect_Spnanob32_ATAC@data[final_list,])

    df.plot$h3k27me3_scale <- (df.plot$h3k27me3 - min(df.plot$h3k27me3)) / (max(df.plot$h3k27me3) - min(df.plot$h3k27me3))
    df.plot$h3k4me3_scale <- (df.plot$h3k4me3 - min(df.plot$h3k4me3)) / (max(df.plot$h3k4me3) - min(df.plot$h3k4me3))
    df.plot$ATAC_scale <- (df.plot$ATAC - min(df.plot$ATAC)) / (max(df.plot$ATAC) - min(df.plot$ATAC))
    
    df.plot <- df.plot[,-3:-5]
    df.plot <- df.plot[df.plot$pt>7,]
    df_long <- pivot_longer(
      data = df.plot[-1],
      cols = -pt, # Select all columns except pseudotime to pivot
      names_to = "condition", # This will be the new column for the condition names
      values_to = "score" # This will be the new column for the scores
    )
    df_long_2 <- df_long[df_long$pt!=0,]
    
    p <- ggplot(data = df_long_2, aes(x = pt, y = score, color = condition)) +
      # geom_point(alpha = 0.5) + # Plots the points with some transparency
      geom_smooth(method = "loess", se = TRUE, aes(fill = condition), alpha = 0.1) + # Adds a smooth line with a confidence interval
      # scale_color_manual(values = c("rna" = "red", "h3k27me3" = "green", "h3k27mac" = "blue")) +
      labs(x = "Pseudo-time", y = "Score", title = "Chromatin opening") +
      theme_minimal() + # Sets a minimal theme
      theme(legend.title = element_blank())+ # Removes the legend title
      labs(title=gene)
    
    print(p)
    
    pdf(paste0("./new_combine_plot/", outputname,"_lineage_noRNA_normalize_", gene, ".pdf"), width = 6.5, height = 6)
    print(p)
    dev.off()
  }
}

plot_one_gene_lineage_RNA_nor <- function(gene){
  final_list <- gene

  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$rna_log@data), final_list)
  
  if (length(final_list) !=0){
    df.plot <- data.frame('cell' = rownames(seurat.wnn.peak.sub@meta.data),
                          'pt' = seurat.wnn.peak.sub$pt,
                          'rna'= seurat.wnn.peak.sub@assays$rna_log@data[final_list,])
    df.plot$RNA_scale <- (df.plot$rna - min(df.plot$rna)) / (max(df.plot$rna) - min(df.plot$rna))
    
    df.plot <- df.plot[,-3]
    
    df.plot <- df.plot[df.plot$pt>7,]
    df_long <- pivot_longer(
      data = df.plot[-1],
      cols = -pt, # Select all columns except pseudotime to pivot
      names_to = "condition", # This will be the new column for the condition names
      values_to = "score" # This will be the new column for the scores
    )
    
    df_long_2 <- df_long[df_long$pt!=0,]
    
    p <- ggplot(data = df_long_2, aes(x = pt, y = score, color = condition)) +
      # geom_point(alpha = 0.5) + # Plots the points with some transparency
      geom_smooth(method = "loess", se = TRUE, aes(fill = condition), alpha = 0.1) + # Adds a smooth line with a confidence interval
      # scale_color_manual(values = c("rna" = "red", "h3k27me3" = "green", "h3k27mac" = "blue")) +
      labs(x = "Pseudo-time", y = "Score", title = "Chromatin opening") +
      theme_minimal() + # Sets a minimal theme
      theme(legend.title = element_blank())+ # Removes the legend title
      labs(title=gene)
    
    print(p)
    
    pdf(paste0("./new_combine_plot/", outputname,"_lineage_RNA_normalize_", gene, ".pdf"), width = 6.5, height = 6)
    print(p)
    dev.off()
  }
}


gene_to_look <- c('Sox2', 'Car10', 'Ank3', 'Gria2')
for (gene_name in gene_to_look){
  plot_one_gene_lineage_no_RNA_nor(c(gene_name))
}

for (gene_name in gene_to_look){
  plot_one_gene_lineage_RNA_nor(c(gene_name))
}



# bivalent (this is for nano68, fig3, you should edit the code for each sample)
plot_one_gene_lineage_all_bi_nor <- function(seurat.wnn.peak.sub, gene){
  final_list <- gene
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore_intersect_Spananob68_deep_H3K4me3@data), final_list)
  final_list <- intersect(rownames(seurat.wnn.peak.sub@assays$gscore_intersect_Spananob68_deep_H3K27me3@data), final_list)
  
  if (length(final_list) !=0){
    df.plot <- data.frame('cell' = rownames(seurat.wnn.peak.sub@meta.data),
                          'pt' = seurat.wnn.peak.sub$Lineage1,
                          'H3K4me3' = seurat.wnn.peak.sub@assays$gscore_intersect_Spananob68_deep_H3K4me3@data[final_list,],
                          'H3K27me3' = seurat.wnn.peak.sub@assays$gscore_intersect_Spananob68_deep_H3K27me3@data[final_list,], 
                          'RNA' = seurat.wnn.peak.sub@assays$SCT@data[final_list,])
    # 'rna'= seurat.wnn.peak.sub@assays$rna_log@data[final_list,])
    df.plot$H3K4me3_scale <- (df.plot$H3K4me3 - min(df.plot$H3K4me3)) / (max(df.plot$H3K4me3) - min(df.plot$H3K4me3))
    df.plot$H3K27me3_scale <- (df.plot$H3K27me3 - min(df.plot$H3K27me3)) / (max(df.plot$H3K27me3) - min(df.plot$H3K27me3))
    df.plot$RNA_scale <- (df.plot$RNA - min(df.plot$RNA)) / (max(df.plot$RNA) - min(df.plot$RNA))
    
    df.plot$bivalency_score <- (df.plot$H3K4me3 + df.plot$H3K27me3) / (abs(df.plot$H3K4me3 - df.plot$H3K27me3) + 1)
    df.plot$bivalency_score_nor <- (df.plot$bivalency_score - min(df.plot$bivalency_score)) / (max(df.plot$bivalency_score) - min(df.plot$bivalency_score))
    
    df.plot <- df.plot[, c('cell', 'pt', 'bivalency_score_nor', 'H3K4me3_scale', 'H3K27me3_scale', 'RNA_scale')]
    
    
    df_long <- pivot_longer(
      data = df.plot[-1],
      cols = -pt, # Select all columns except pseudotime to pivot
      names_to = "condition", # This will be the new column for the condition names
      values_to = "score" # This will be the new column for the scores
    )
    
    df_long_2 <- df_long[df_long$pt!=0,]
    
    p <- ggplot(data = df_long_2, aes(x = pt, y = score, color = condition)) +
      # geom_point(alpha = 0.5) + # Plots the points with some transparency
      geom_smooth(method = "loess", se = TRUE, aes(fill = condition), alpha = 0.1) + # Adds a smooth line with a confidence interval
      # scale_color_manual(values = c("rna" = "red", "h3k27me3" = "green", "h3k27mac" = "blue")) +
      labs(x = "Pseudo-time", y = "Score", title = "Chromatin opening") +
      theme_minimal() + # Sets a minimal theme
      theme(legend.title = element_blank())+ # Removes the legend title
      labs(title=gene)
    
    print(p)
    
    pdf(paste0("./new_combine_plot/", outputname,"_lineage_all_bi-no-nor_normalize_", gene, ".pdf"), width = 6.5, height = 6)
    print(p)
    dev.off()
  }
}


plot_one_gene_lineage_all_bi_nor(seurat.wnn.peak.sub.clean, 'Sox2')


