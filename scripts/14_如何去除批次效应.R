### 批次效应的本质是不同批次数据因技术差异（如测序平台、实验时间、操作人员等）导致的系统性偏差，这些偏差可能掩盖真实的生物差异（如细胞类型、状态的差异）。rPCA通过正则化项约束批次相关的变异，减少批次效应对主成分（PCs）的影响，同时保留生物相关的变异。其核心是：对不同批次数据的协方差结构引入正则化调整，降低批次特异性变异的权重；使校正后的主成分更能反映真实的生物异质性，而非技术批次差异。

library(data.table)
library(tidyverse)
library(ggpubr)
library(Seurat)
scRNAa<-Read10X("scRNA_read/500_PBMC_3p_LT_Chromium_X_filtered_feature_bc_matrix/filtered_feature_bc_matrix/")


#### MTX
library(Matrix)
# Read in `matrix.mtx`
counts <- readMM("scRNA_read/frozen_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/matrix.mtx")
dim(counts)
counts[1:3,1:3]
# Read in `genes.tsv`
library(readr)
genes <- read_tsv("scRNA_read/frozen_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/genes.tsv", col_names = FALSE)
gene_ids <- genes$X2
# Read in `barcodes.tsv`
cells <- read_tsv("scRNA_read/frozen_pbmc_donor_a_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/barcodes.tsv", col_names = FALSE)

cell_ids <- cells$X1
# Create a sparse matrix for more efficient computation
counts <- as(counts, "dgCMatrix")
#counts <-as.data.frame(counts)
# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts[1:3,1:3]
head(counts)

colnames(counts)<-paste0("A_",colnames(counts))
rownames(counts)<-gsub("_","-",rownames(counts))
scRNAb=counts

###############################################################################
### 走scRNA流程

sca <- CreateSeuratObject(scRNAa,project = "A")
sca$group="A"
scb <- CreateSeuratObject(scRNAb,project = "B")
scb$group="B"
###############################################################################
########## 不做任何处理
batch_sc <-merge(sca,scb)

batch_sc <- NormalizeData(batch_sc)
batch_sc <- FindVariableFeatures(batch_sc, selection.method = "vst", nfeatures = 2000)

batch_sc <- ScaleData(batch_sc, verbose = FALSE)
batch_sc <- RunPCA(batch_sc, npcs = 30, verbose = FALSE)
batch_sc <- RunUMAP(batch_sc, reduction = "pca", dims = 1:30)
batch_sc <- FindNeighbors(batch_sc, reduction = "pca", dims = 1:30)
batch_sc <- FindClusters(batch_sc, resolution = 0.5)

DimPlot(batch_sc,group.by = "group")


###############################################################################
########## Seurat 推荐 rPCA（Regularized PCA，正则化主成分分析），去除批次效应的统计方法。

sc.list <- list(a=sca,
                b=scb)

# bize and identify variable features for each dataset independently
sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = sc.list)

sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = sc.list, 
                                  anchor.features = features, 
                                  reduction = "rpca")

# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

DimPlot(combined,group.by = "group")


####################################################################################
##### harmony
library(harmony)
harmony_sc <-merge(sca,scb)

harmony_sc <- NormalizeData(harmony_sc)
harmony_sc <- FindVariableFeatures(harmony_sc, selection.method = "vst", nfeatures = 2000)

harmony_sc <- ScaleData(harmony_sc, verbose = FALSE)
harmony_sc <- RunPCA(harmony_sc, npcs = 30, verbose = FALSE)
harmony_sc$group
harmony_sc <- harmony_sc %>% 
  RunHarmony("group", plot_convergence = TRUE)

harmony_sc <- harmony_sc %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(harmony_sc,group.by = "group")

################################################################################