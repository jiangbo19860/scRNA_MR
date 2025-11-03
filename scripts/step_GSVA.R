rm(list = ls())

library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(GSVA)

gmtfile <- "h.all.v7.4.symbols.gmt"

hallmark <- read.gmt(gmtfile)

hallmark$term <- gsub('HALLMARK_','',hallmark$term)

hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)


load("scRNA_demo.rds")
DimPlot(scRNA)
expr <- as.matrix(scRNA@assays$RNA@counts)

dim(expr)

es.matrix = gsva(expr, 
                 hallmark.list, 
                 kcdf="Poisson",
                 method="ssgsea", 
                 abs.ranking=T ,
                 parallel.sz=3)


es.matrix =as.data.frame(t(es.matrix))

scRNA <-AddMetaData(scRNA,es.matrix)

meta <-scRNA@meta.data

colnames(meta)

FeaturePlot(scRNA,features = c("TGF_BETA_SIGNALING"))
