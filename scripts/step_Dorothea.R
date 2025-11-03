BiocManager::install("dorothea")
library(tidyverse)
library(pheatmap)
## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

library(dorothea)
library(Seurat)
load("scRNA_endo.Rds")
scRNA <- run_viper(scRNA, regulon,
                  options = list(method = "scale", 
                                 minsize = 4, 
                                 eset.filter = FALSE, 
                                 cores = 1, 
                                 verbose = FALSE))

DimPlot(scRNA,group.by = "RNA_snn_res.0.1")
## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = scRNA) <- "dorothea"
scRNA <- ScaleData(scRNA)

scRNA.TF<-FindAllMarkers(scRNA,only.pos = T)
scRNA.TF.top<- scRNA.TF %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
DoHeatmap(scRNA,features = scRNA.TF.top$gene)
# ClusterName
Idents(scRNA)<-"ClusterName"

scRNA.TF<-FindAllMarkers(scRNA,only.pos = T)
scRNA.TF.top<- scRNA.TF %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
DoHeatmap(scRNA,features = scRNA.TF.top$gene)


# CellFromTumor

Idents(scRNA)<-"CellFromTumor"

scRNA.TF<-FindAllMarkers(scRNA,only.pos = T)
scRNA.TF.top<- scRNA.TF %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
DoHeatmap(scRNA,features = scRNA.TF.top$gene)



















