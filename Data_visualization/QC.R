####### This code is for QC analysis

library(ggplot2)
library(dplyr)

###### ATAC + Histone modification
## Unique fragments
# Read the data from the TXT files
setwd("/mnt/nas1/Users/Pengfei/Processed_data/2023/Nanobody_project/meta_data")
files <- list.files()
data_list <- lapply(files, function(file) {
  print(file)
  read.table(file, header = TRUE, sep = " ")
})

data_list[[1]]
length(data_list)

# Combine the data frames and add a 'Source' column
new_files <- gsub("_meta_data.txt", "", files)

for (i in 1:length(data_list)){
  file_name <- new_files[i]
  data_list[[i]]$Source <- file_name
}

# compare 50um
setwd("/mnt/nas1/Users/Liran/Processed_data/2024/QC")

combined_data <- dplyr::bind_rows(data_list[1], data_list[2],data_list[3], data_list[4],data_list[5], data_list[6],data_list[7])
colnames(combined_data)
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 7, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(1.22, 0.5), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"),
          axis.text.x = element_text(angle = 20, hjust = 1) 
    )
}

distributions2 <- ggplot(combined_data, aes(x = Source, y = log10(nFrags))) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()
distributions2

pdf(paste0("QC_paper/50um_qc_log10nFrags.pdf"), width = 9, height = 5)
print(distributions2)
dev.off()


distributions2 <- ggplot(combined_data, aes(x = Source, y =TSSEnrichment)) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()

distributions2
pdf(paste0("QC_paper/50um_qc_TSSEnrichment.pdf"), width = 9, height = 5)
print(distributions2)
dev.off()

##### compare 20
combined_data <- dplyr::bind_rows(data_list[8], data_list[9],data_list[10], data_list[11])
colnames(combined_data)

distributions2 <- ggplot(combined_data, aes(x = Source, y = log10(nFrags))) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()
distributions2

theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 6, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(1.25, 0.5), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"),
          axis.text.x = element_text(angle = 20, hjust = 1) 
    )
}

pdf(paste0("QC_paper/20um_qc_log10nFrags.pdf"), width = 7, height = 5)
print(distributions2)
dev.off()


distributions2 <- ggplot(combined_data, aes(x = Source, y =TSSEnrichment)) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()

distributions2
pdf(paste0("QC_paper/20um_qc_TSSEnrichment.pdf"), width = 7, height = 5)
print(distributions2)
dev.off()


## fragmentSizes
setwd('/mnt/nas1/Users/Liran/Processed_data/2024/')
samplepath_list <- c('/mnt/nas1/Users/Liran/Processed_data/2024/17nanobody_deep/Spnanob17_h3k27ac_deep/Save-inTissue-Spananob17_H3K27ac',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/17nanobody_deep/Spnanob17_h3k27me3_deep/Save-inTissue-Spananob17_H3K27me3',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/30nanobody_replicate_of_17/Spnanob30_H3K27me3/Save-inTissue-Spananob30_h3k27me3',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/30nanobody_replicate_of_17/Spnanob30_H3K27ac/Save-inTissue-Spananob30_h3k27ac',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/32nanobody/Spnanob32_H3K27me3/Save-inTissue-Spnanob32_H3K27me3',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/32nanobody/Spnanob32_H3K4me3/Save-inTissue-Spnanob32_H3K4me3',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/32nanobody/Spnanob32_ATAC/Save-inTissue-Spnanob32_ATAC',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/38nanobody_deep/Spnanob38_deep_H3K27ac/Save-inTissue-Spananob38_deep_H3K27ac',
                     '/mnt/nas1/Users/Liran/Processed_data/2024/38nanobody_deep/Spnanob38_deep_H3K27me3/Save-inTissue-Spananob38_deep_H3K27me3',
                     '/mnt/nas1/Users/Pengfei/Processed_data/2023/Nanobody_project/Spnb68/Spnanob68_H3K4me3/Spnanob68/fragments/Save-inTissue-Spnanob68_H3K4me3',
                     '/mnt/nas1/Users/Pengfei/Processed_data/2023/Nanobody_project/Spnb68/Spnanob68_H3K27me3/fragments/Save-inTissue-Spnanob68_H3K27me3')


i <- 1
samplepath <- samplepath_list[i]
projCUTA <- loadArchRProject(path = samplepath, force = FALSE, showLogo = TRUE)
p <- plotFragmentSizes(ArchRProj = projCUTA, returnDF = TRUE)
fragmentSizes <- p
for (i in 2:9){
  samplepath <- samplepath_list[i]
  projCUTA <- loadArchRProject(path = samplepath, force = FALSE, showLogo = TRUE)
  fragmentSizes_add <- plotFragmentSizes(ArchRProj = projCUTA, returnDF = TRUE)
  fragmentSizes <- rbind(fragmentSizes, fragmentSizes_add)
  
}


fragmentSizes_df <- data.frame(fragmentSizes)
sample_name <- unique(fragmentSizes_df$group)
plot <- ggplot(fragmentSizes_df, aes(x = fragmentSize , y = fragmentPercent, color = group)) +
  geom_line() + 
  theme_minimal() +
  labs(x = "Fragment Size (bp)", y = "Percentage of Fragments", title = "Fragment Sizes Across Projects") +
  scale_color_brewer(palette = "Paired")  

pdf(paste0("QC_paper/Fragment_Size_all.pdf"), width = 8, height = 6)
print(plot)
dev.off()



# plot per chip-size
# 50um
sample_name
filtered_df <- fragmentSizes_df %>%
  filter(group %in% c("Spananob17_H3K27ac", "Spananob17_H3K27me3", "Spananob30_h3k27me3",     
                      "Spananob30_h3k27ac", "Spnanob32_H3K27me3", "Spnanob32_H3K4me3",       
                      "Spnanob32_ATAC" ))


plot50 <- ggplot(filtered_df, aes(x = fragmentSize , y = fragmentPercent, color = group)) +
  geom_line() + 
  theme_minimal() +
  labs(x = "Fragment Size (bp)", y = "Percentage of fragments", title = "Fragment Sizes Across Projects (50um)") +
  scale_color_brewer(palette = "Set1")  

pdf(paste0("QC_paper/Fragment_Size_50um.pdf"), width = 8, height = 6)
print(plot50)
dev.off()


# 20um
sample_name
filtered_df <- fragmentSizes_df %>%
  filter(group %in% c("Spnanob68_H3K4me3", "Spnanob68_H3K27me3", "Spananob38_deep_H3K27ac", "Spananob38_deep_H3K27me3" ))


plot20 <- ggplot(filtered_df, aes(x = fragmentSize , y = fragmentPercent, color = group)) +
  geom_line() + 
  theme_minimal() +
  labs(x = "Fragment Size (bp)", y = "Percentage of fragments", title = "Fragment Sizes Across Projects (20um)") +
  scale_color_brewer(palette = "Set1")  # Adjust colors as necessary

pdf(paste0("QC_paper/Fragment_Size_20um.pdf"), width = 8, height = 6)
print(plot20)
dev.off()


###### RNA
# Read the data from the TXT files
setwd("/mnt/nas1/Users/Liran/Processed_data/2024/QC/RNA")
files <- list.files()
spatial.obj@meta.data

meta_save <- list()

for (file in files){
  spatial <- readRDS(file)
  meta_save[[file]]<-spatial@meta.data
}

for (file in files){
  meta_save[[file]]$Source <- file
}

# compare 50um
combined_data <- dplyr::bind_rows(meta_save[['32_seurat.rna.rds']], meta_save[['38_seurat.rna.rds']],meta_save[["68_seurat.rna.rds"]])

colnames(combined_data)
dim(combined_data)
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 7, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(1.22, 0.5), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"),
          axis.text.x = element_text(angle = 20, hjust = 1) 
    )
}

distributions2 <- ggplot(combined_data, aes(x = Source, y = log10(nCount_Spatial))) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()
distributions2

pdf(paste0("../QC_paper/RNA_qc_log10nFrags.pdf"), width = 8, height = 5)
print(distributions2)
dev.off()


distributions2 <- ggplot(combined_data, aes(x = Source, y =nFeature_Spatial)) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()

distributions2
pdf(paste0("../QC_paper/RNA_qc_nFeature.pdf"), width = 8, height = 5)
print(distributions2)
dev.off()


distributions2 <- ggplot(combined_data, aes(x = Source, y =log10(nFeature_Spatial))) +
  geom_violin(aes(fill = Source, colour = Source), alpha = 0.5) +
  geom_boxplot(aes(colour = Source), width = 0.2) +
  theme_niwot()

distributions2
pdf(paste0("../QC_paper/RNA_qc_log10nFeature.pdf"), width = 8, height = 5)
print(distributions2)
dev.off()
