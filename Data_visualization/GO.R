####### This code is for GO analysis

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# markerGenes_go is the list pf genes for GO analysis
gene.100<-data.frame("SYMBOL"=toupper(markerGenes_go))
gene.100.enterzid <- bitr(gene.100$SYMBOL,fromType = 'SYMBOL',toType = c('ENSEMBL'),OrgDb = org.Hs.eg.db)
gene.100.go<-enrichGO(gene = gene.100.enterzid$ENSEMBL,OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",ont="ALL", pAdjustMethod="BH",pvalueCutoff = 0.05, qvalueCutoff = 0.2) 

# plot and save
go_dot<-dotplot(gene.100.go, split="ONTOLOGY",font.size=10,label_format = 50)+facet_grid(ONTOLOGY~., scale="free")
pdf("./new_combine_plot/GO_C8.pdf",width=10,height=10)
print(go_dot)
dev.off()

