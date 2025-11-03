# install.packages("devtools")
  devtools::install_github("saeyslab/nichenetr")
  library(nichenetr)
  library(Seurat)
  library(tidyverse)
  rm(list=ls())
  
  ## 读取NicheNet先验数据
  ligand_target_matrix <- readRDS("nichenetr/ligand_target_matrix.rds")
  lr_network <- readRDS("nichenetr/lr_network.rds")
  weighted_networks <- readRDS("nichenetr/weighted_networks.rds")
  
  ## 读取单细胞数据
  scRNA <-Read10X_h5("TISCH/BLCA_GSE130001_expression.h5")
  scRNA <-CreateSeuratObject(counts = scRNA)
  meta <-data.table::fread("TISCH/BLCA_GSE130001_CellMetainfo_table.tsv",
                           header = T,
                           data.table = F)
  rownames(meta)<-meta$Cell
  scRNA<-AddMetaData(scRNA,metadata = meta)
  
  scRNA <- NormalizeData(scRNA, 
                         normalization.method = "LogNormalize", 
                         scale.factor = 10000) 
  
  
  ##选择高变基因
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(scRNA), 10) 
  plot1 <- VariableFeaturePlot(scRNA) 
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
  CombinePlots(plots = list(plot1, plot2),legend="bottom") 
  
  ##数据中心化与回归
  #数据中心化
  scale.genes <-  rownames(scRNA)
  #scale.genes <-  VariableFeatures(scRNA)    #内存不够可用此命令代替上一行
  scRNA <- ScaleData(scRNA, features = scale.genes)
  ##提取主成分
  #PCA降维去噪
  scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
  
  plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
  
  plot2 <- ElbowPlot(scRNA, ndims=30, reduction="pca") 
  plot1+plot2
  
  #选取主成分
  pc.num=1:20
  
  ##细胞聚类
  scRNA <- FindNeighbors(scRNA, dims = pc.num) 
  scRNA <- FindClusters(scRNA, resolution =0.5)
  table(scRNA@meta.data$seurat_clusters)
  metadata <- scRNA@meta.data
  
  ##非线性降维
  #tSNE
  scRNA <- RunTSNE(scRNA, dims = pc.num)
  DimPlot(scRNA, reduction = "tsne", label=T) 
  DimPlot(scRNA, reduction = "tsne", label=T,group.by = "Celltype..minor.lineage.")
  #UMAP
  scRNA <- RunUMAP(scRNA, dims = pc.num)
  DimPlot(scRNA, reduction = "umap", label=T) 
  
  ## 一步完成单细胞数据的NicheNet分析
  Idents(scRNA) <- "Celltype..minor.lineage."
  unique(scRNA$Celltype..minor.lineage.)
  sender.cells <- c("Endothelial","Fibroblasts","Myofibroblasts")
  table(metadata$Patient)
  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = scRNA, 
                                                 top_n_ligands = 20,
                                                 receiver = "Epithelial", 
                                                 sender = sender.cells,
                                                 condition_colname = "Patient", 
                                                 condition_oi = "sample1", 
                                                 condition_reference = "sample2", 
                                                 ligand_target_matrix = ligand_target_matrix, 
                                                 lr_network = lr_network, 
                                                 weighted_networks = weighted_networks, 
                                                 organism = "human")
  
  

  scRNA$celltype_group <-paste0(scRNA$Celltype..minor.lineage., ": ", scRNA$Patient) 
  
  ### 查看top配体在各个细胞亚型的表达情况
  p = DotPlot(subset(scRNA, Celltype..minor.lineage. %in% sender.cells), 
              group.by = "celltype_group", 
              features = nichenet_output$top_ligands, 
              cols = "RdYlBu") + RotatedAxis()
  p
  ### 查看top受体在Epithelial细胞中的表达情况
  p = VlnPlot(subset(scRNA, Celltype..minor.lineage.=="Epithelial"), 
              features = nichenet_output$top_receptors[1:10], 
              group.by = "Patient", 
              split.by = "Patient", ncol = 5)
  p
  
  ### 查看top靶基因在Epithelial细胞中的表达情况
  p = VlnPlot(subset(scRNA, Celltype..minor.lineage.=="Epithelial"), 
              features = nichenet_output$top_targets[1:10], 
              group.by = "Patient", 
              split.by = "Patient", ncol = 5, 
              combine = T)
  p
  
  
  ### 配体与靶基因的调控关系
  p <- nichenet_output$ligand_target_heatmap + theme(legend.position = "right")
  
  p <- nichenet_output$ligand_activity_target_heatmap
  
  df <-nichenet_output$ligand_target_df
  
  write.csv(df, "ligand_target_df.csv", row.names = F)
  
  ### 配体与受体的互作关系
  p <- nichenet_output$ligand_receptor_heatmap
  
  df <- nichenet_output$ligand_receptor_df
  
  write.csv(df, "ligand_receptor_df.csv", row.names = F)
  
  # 'bona fide'是指有文献报道的配体-受体，而不是基于蛋白质互作（PPI）预测的。
  p <- nichenet_output$ligand_receptor_heatmap_bonafide
  ggsave("ligand_receptor_heatmap_bonafide.pdf", p, width = 12, height = 6.5)
  
  df <- nichenet_output$ligand_receptor_df_bonafide
  write.csv(df, "ligand_receptor_df_bonafide.csv", row.names = F)
