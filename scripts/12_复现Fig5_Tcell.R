rm(list = ls())
pacman::p_load(here, Seurat, dplyr, data.table, ggplot2, cowplot, patchwork, ggrepel, hdf5r, Matrix, ggpubr, forcats, BiocParallel, pheatmap, RColorBrewer, GSVA, clusterProfiler, limma, stringr, pheatmap)

here()

load(here("1_data/E-MTAB-6149/T_cell.Cellview.Rds"))
log2cpm[1:3, 1:3]
load(here("1_data/E-MTAB-6149/scRNA.Rds"))
dim(scRNA) # 行gene，列cell

scRNA <- scRNA[, colnames(log2cpm)]
scRNA <- AddMetaData(scRNA, metadata = tsne.data)
# 使用显式的 layer 参数（适用于 Seurat 5+）
scRNA <- NormalizeData(scRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  assay = DefaultAssay(scRNA), # 指定要归一化的分析
  layer = "counts"
)

## 从表达矩阵中识别高变基因（Highly Variable Genes, HVGs）
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(scRNA), 10)
## 数据中心化与回归
# 数据中心化
# scale.genes <-  rownames(scRNA)
# 获取Seurat对象中已识别的高变基因列表
scale.genes <- VariableFeatures(scRNA)

scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

ElbowPlot(scRNA, ndims = 20, reduction = "pca")

pc.num <- 1:15 # 拐点出现在 PC15，表示前 15 个 PC 捕获了大部分方差）。

scRNA <- FindNeighbors(scRNA, dims = pc.num)

scRNA <- FindClusters(scRNA, resolution = 0.2)

metadata <- scRNA@meta.data
colnames(metadata)

scRNA <- RunTSNE(scRNA, dims = pc.num)

scRNA <- RunUMAP(scRNA, dims = pc.num) # R 原生的 uwot 包，使用余弦距离（cosine）作为度量标准。

save(scRNA, file = "scRNA_T_cell.Rds")

load("scRNA_T_cell.Rds")

# 打印 scRNA 元数据的列名
colnames(scRNA@meta.data)

DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne")

DimPlot(scRNA, group.by = "ClusterName", reduction = "tsne")

DimPlot(scRNA, reduction = "tsne", label = TRUE) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.position = "none"
  ) +
  labs(title = "T cell") +
  ggsci::scale_color_lancet()



ggsave(filename = "Fig/Fig5a1按seurat_clusters分簇的t-SNE图.pdf", height = 8, width = 8)

dim(scRNA)
unique(scRNA$seurat_clusters)

################################################################
# Fig5a2_按细胞来源CellFromTumor分组的t-SNE图 ------------------------------------
DimPlot(scRNA,
  reduction = "tsne",
  group.by = "CellFromTumor"
) +
  scale_color_manual(
    values = c("#72C667", "#2D5474"),
    labels = c("Non-malignant", "Tumor")
  ) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    legend.position = c(.01, .1)
  ) +
  labs(title = "Sample Origin")

ggsave(filename = "Fig/Fig5a2_按细胞来源CellFromTumor分组的t-SNE图.pdf", height = 8, width = 8)

################################################################
# Fig1d_T cell水平堆叠柱状图 ------------------------------------------------------
tmp <- scRNA@meta.data

# 1. 处理 seurat_clusters 因子水平
sorted_clusters <- tmp$seurat_clusters %>%
  as.character() %>%
  as.numeric() %>%
  unique() %>%
  sort(decreasing = TRUE) %>%
  as.character()
tmp$seurat_clusters <- factor(tmp$seurat_clusters, levels = sorted_clusters)

# 2. group 因子水平：levels = c("Tumor", "Non-malignant")
tmp$group <- factor(
  ifelse(tmp$CellFromTumor == "1", "Tumor", "Non-malignant"),
  levels = c("Tumor", "Non-malignant")
)

# 3. 绘图
p1 <- ggplot(tmp, aes(x = seurat_clusters, fill = group)) +
  geom_bar(stat = "count", position = "fill") + # 按因子水平顺序堆叠
  scale_fill_manual(
    values = c("#2D5474", "#72C667"),
    labels = c("Tumor", "Non-malignant")
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = "black") # 可选：显示 y 轴聚类标签
  ) +
  xlab("") +
  ylab("") +
  coord_flip()

p1

# p2 绘图代码同理
p2 <- ggplot(
  tmp,
  aes(x = seurat_clusters, fill = PatientNumber)
) +
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  scale_fill_manual(values = c("#268A24", "#8EEE8B", "#FB6346", "#FBD51A", "#28507D")) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "top",
    axis.line.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab(" ") +
  ylab(" ") +
  coord_flip()

# 组合绘图
p1 + p2 + plot_layout(ncol = 2)
ggsave("Fig/Fig1d_T_cell_水平堆积柱状图.pdf", height = 4, width = 10)

################################################################
# Fig5b_T cell特征基因的t-SNE图 ---------------------------------------------------
p <- FeaturePlot(scRNA,
  reduction = "tsne",
  features = c(
    "FGFBP2",
    "CD8A",
    "CD4",
    "FOXP3",
    "LAG3",
    "MKI67"
  ),
  ncol = 3,
  cols = c("gray", "red")
)
p & theme(
  panel.border = element_rect(
    fill = NA, color = "black",
    size = 1, linetype = "solid"
  ),
  legend.position = "none"
)

ggsave("Fig/Fig5b_Tcell特征基因的t-SNE图.pdf", height = 6, width = 9)

################################################################
# Fig5d_热图pheatmap展示用GSVA（基因集变异分析） 和limma差异分析不同T细胞亚群的特定通路在肿瘤细胞vs非肿瘤细胞中的通路活性差异-----------------
Idents(scRNA) <- "ClusterName"
DimPlot(scRNA, reduction = "tsne", label = T)

library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

gs <- read.gmt(here("1_data/GSEA/h.all.v7.4.symbols.gmt"))

gs$term <- gsub("HALLMARK_", "", gs$term)

gs.list <- gs %>%
  split(.$term) %>%
  lapply("[[", 2)

# 分群并循环处理每个细胞亚群
clusternames <- unique(scRNA@active.ident) # 获取所有细胞亚群名称

data_tmp <- data.frame(row.names = unique(gs$term))

for (i in 1:4) {
  scRNAsub <- subset(scRNA, ClusterName == clusternames[i]) # 提取特定亚群的细胞
  expr_norm <- GetAssayData(scRNAsub, slot = "data") # 替代 counts，使用归一化后的数据
  register(MulticoreParam(workers = 3)) # 3个核心，可根据电脑配置调整
  ssgsea_param <- ssgseaParam(
    expr = expr_norm, # 输入归一化后的表达矩阵（基因×细胞）
    geneSets = gs.list # 输入基因集列表
  )
  gsva_es <- gsva(
    param = ssgsea_param, # 传入参数对象
    BPPARAM = bpparam(), # 若已设置并行计算，会自动应用（如之前的 register(MulticoreParam)）
    verbose = FALSE # 可选：静默运行，不输出中间信息
  )
  # 分组（肿瘤 vs 非肿瘤）
  scRNAsub$group <- ifelse(scRNAsub$CellFromTumor == "1", "tumor", "Non_malignant")

  group_list <- data.frame(
    sample = colnames(gsva_es),
    group = scRNAsub$group
  )

  head(group_list)
  # 差异分析（limma 包）构建实验设计矩阵（Design Matrix）目的是将分组信息（肿瘤 vs 非肿瘤）转换为 数学模型
  design <- model.matrix(~ 0 + factor(group_list$group)) # ~ 表示 “取决于”“关于” 的关系，分隔模型中的响应变量（~左边）和解释变量，0表示不包含截距项（Intercept），强制每个分组单独作为一个变量。效果：直接比较组间差异（如 tumor - Non_malignant），避免截距项的干扰。生成一个设计矩阵，每行对应一个样本，每列对应一个分组。矩阵中的 1 表示样本属于该分组，0 表示不属于。核心是用 model.matrix() 生成“仅包含解释变量” 的矩阵（用于后续的差异分析建模），因此不需要指定响应变量。
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  design
  contrast.matrix <- makeContrasts(tumor - Non_malignant, levels = design)
  fit <- lmFit(gsva_es, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2,
    coef = 1, n = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  x <- data.frame(
    row.names = rownames(x),
    t = x$t
  )
  colnames(x) <- clusternames[i]
  data_tmp <- cbind(data_tmp, x)
}


write.csv(data_tmp, "T_cell_gsva_limma.csv", quote = F)

head(data_tmp)
str(data_tmp)

# 确保行名（通路名称）和列名（细胞亚群名称）正确
# 可以根据实际情况调整列名的显示，比如将列名中的下划线等替换成空格等更美观的形式
colnames(data_tmp) <- c("CD4+ T cells", "CD8+ T cells", "Regulatory T cells", "Natural killer cells")

# 转换为矩阵（pheatmap 需要矩阵输入）
data_matrix <- as.matrix(data_tmp)

# 检查数据矩阵的基本属性
print(dim(data_matrix)) # 应输出行和列数（如 50 行 4 列）
print(any(is.na(data_matrix))) # 应输出 FALSE（无缺失值）
print(class(data_matrix)) # 应输出 "matrix" 或 "numeric"

# 定义颜色梯度
color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50) # 控制颜色梯度的细腻程度。数值越大，颜色过渡越平滑

# 移除手动打开/关闭PDF设备的代码，直接在pheatmap中指定输出文件
heatmap_params <- list(
  mat = data_matrix,
  color = color_palette,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Pathway activity (t value tumour versus non-malignant)",
  fontsize_row = 8,
  fontsize_col = 10,
  legend = TRUE,
  border_color = NA,
  cellwidth = 30,
  cellheight = 12,
  angle_col = "45",
  filename = here("Fig/Fig5d_pheatmap.pdf"),  # 直接指定输出路径
  height = 12,  # 图片高度（英寸）
  width = 10    # 图片宽度（英寸）
)

# 直接绘图并保存，无需手动调用pdf()和dev.off()
do.call(pheatmap, heatmap_params)

##### CD8+ T cells ###########################################################
# Fig5e_用Seurat、GSVA、limma和pheatmap等工具，完成CD8+ T cells单细胞数据处理 → 分群鉴定 → 通路活性分析

meta <- scRNA@meta.data
unique(meta$ClusterName)

scRNA <- subset(scRNA, ClusterName == "CD8+ T cells ")

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
## 选择高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# top10 <- head(VariableFeatures(scRNA), 10)
## 数据中心化与回归
# 数据中心化
# scale.genes <-  rownames(scRNA)
scale.genes <-  VariableFeatures(scRNA)    #内存不够可用此命令代替上一行
scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

ElbowPlot(scRNA, ndims = 20, reduction = "pca")

pc.num <- 1:15

scRNA <- FindNeighbors(scRNA, dims = pc.num)

scRNA <- FindClusters(scRNA, resolution = 0.1)

metadata <- scRNA@meta.data

scRNA <- RunTSNE(scRNA, dims = pc.num)

scRNA <- RunUMAP(scRNA, dims = pc.num)

save(scRNA, file = "scRNA_CD8T_cell.Rds")

DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne") # cluster 2,4,5,8

DimPlot(scRNA, reduction = "tsne")

# library(GSVA)
# library(clusterProfiler)
# library(limma)
# library(stringr)
# library(ggplot2)

gs <- read.gmt(here("1_data/GSEA/h.all.v7.4.symbols.gmt"))

gs$term <- gsub("HALLMARK_", "", gs$term)

gs.list <- gs %>%
  split(.$term) %>%
  lapply("[[", 2)

# 分群并循环处理每个细胞亚群
clusternames <- unique(scRNA$dbCluster) # 获取所有细胞亚群名称,2,8,5,4

data_tmp <- data.frame(row.names = unique(gs$term))

for (cluster in clusternames) {
  scRNAsub <- subset(scRNA, dbCluster == cluster) # 提取特定亚群的细胞
  # 打印当前处理的亚群
  message("===== Processing cluster: ", cluster, " =====")
  # 检查分组是否完整（必须同时有tumor和Non_malignant）
  group_counts <- table(scRNAsub$CellFromTumor)
  message("  Group counts (1=tumor, 0=Non_malignant): ")
  print(group_counts)
  expr_norm <- GetAssayData(scRNAsub, slot = "data") # 替代 counts，使用归一化后的数据
  register(MulticoreParam(workers = 3)) # 3个核心，可根据电脑配置调整
  ssgsea_param <- ssgseaParam(
    expr = expr_norm, # 输入归一化后的表达矩阵（基因×细胞）
    geneSets = gs.list # 输入基因集列表
  )
  gsva_es <- gsva(
    param = ssgsea_param, # 传入参数对象
    BPPARAM = bpparam(), # 若已设置并行计算，会自动应用（如之前的 register(MulticoreParam)）
    verbose = FALSE # 可选：静默运行，不输出中间信息
  )
  # 分组（肿瘤 vs 非肿瘤）
  scRNAsub$group <- ifelse(scRNAsub$CellFromTumor == "1", "tumor", "Non_malignant")
  
  group_list <- data.frame(
    sample = colnames(gsva_es),
    group = scRNAsub$group
  )
  
  head(group_list)
  
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  design
  contrast.matrix <- makeContrasts(tumor - Non_malignant, levels = design)
  fit <- lmFit(gsva_es, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2,
    coef = 1, n = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  x <- data.frame(
    row.names = rownames(x),
    t = x$t
  )
  colnames(x) <- paste0("CLuster", cluster)
  data_tmp <- cbind(data_tmp, x)
}

# 检查每个亚群的细胞数量
table(scRNA$dbCluster)

# 可视化某个亚群的基因表达分布
VlnPlot(scRNA, features = "CD3E", group.by = "dbCluster")


# 转换为矩阵（pheatmap 需要矩阵输入）
data_matrix <- as.matrix(data_tmp)

# 检查数据矩阵的基本属性
print(dim(data_matrix)) # 应输出行和列数（如 50 行 4 列）
print(any(is.na(data_matrix))) # 应输出 FALSE（无缺失值）
print(class(data_matrix)) # 应输出 "matrix" 或 "numeric"

# 定义颜色梯度
color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50) # 控制颜色梯度的细腻程度。数值越大，颜色过渡越平滑

# 移除手动打开/关闭PDF设备的代码，直接在pheatmap中指定输出文件
heatmap_params <- list(
  mat = data_matrix,
  color = color_palette,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Pathway activity (t value tumour versus non-malignant)",
  fontsize_row = 8,
  fontsize_col = 10,
  legend = TRUE,
  border_color = NA,
  cellwidth = 30,
  cellheight = 12,
  angle_col = "45",
  filename = here("Fig/Fig5d_CD8T_pheatmap.pdf"),  # 直接指定输出路径
  height = 12,  # 图片高度（英寸）
  width = 10    # 图片宽度（英寸）
)

# 直接绘图并保存，无需手动调用pdf()和dev.off()
do.call(pheatmap, heatmap_params)


################################################################################
# 堆叠小提琴图 ------------------------------------------------------------------
gene_name <- c(
  "PDCD1",
  "CD28",
  "CTLA4",
  "ICOS",
  "BTLA",
  "LAG3",
  "TNFRSF9",
  "TNFRSF4",
  "CD27",
  "HAVCR2",
  "GZMA",
  "GZMB",
  "GZMM",
  "GZMH",
  "GZMK"
)
source(here("2_src/编程猫_src/stackvlion.R"))
# StackedVlnPlot(scRNA, features = gene_name, idents = c(0, 1, 2, 3)) & ggsci::scale_fill_lancet()

StackedVlnPlot(scRNA, features = gene_name, idents = c(0,1,2,3)) +
  ggsci::scale_fill_lancet()  +
  scale_y_continuous(limits = c(0, 5))  # 手动指定Y轴范围（根据实际数据调整）

ggsave(filename = here("Fig/Fig5e2_堆积小提琴图.pdf"), height = 15, width = 6)


################################################################################

# 颗粒酶（GZMA/GZMB/GZMH）相关基因的全基因组相关性分析与可视化。基因表达可视化→数据过滤→综合评分构建→全基因组 Spearman 相关性分析→结果可视化。系统筛选与颗粒酶（GZMA、GZMB、GZMH）表达相关的基因，聚焦最显著的前 30 个基因。解析免疫细胞的细胞毒性功能。
# 关闭所有图形设备，重新绘图
while (dev.cur() > 1) dev.off()
FeaturePlot(
  object = scRNA,
  features = c("GZMA", "GZMB", "GZMH"),
  reduction = "umap"  # 明确指定已存在的降维方法
)


exp <- GetAssayData(scRNA)
exp <- as.data.frame(exp)
exp <- exp[, colSums(exp) > 0]
exp <- exp[rowSums(exp) > 0, ]

GZ <- exp[c("GZMA", "GZMB", "GZMH"), ]
GZ <- as.data.frame(GZ)
for (i in 1:ncol(GZ)) {
  GZ["GZ", i] <- mean(as.numeric(GZ[1:3, i]))
}

y <- as.numeric(GZ[4, ])
rownames <- rownames(exp)
cor <- do.call(rbind, lapply(rownames, function(x) {
  dd <- cor.test(as.numeric(exp[x, ]), y, type = "spearman")
  data.frame(
    gene = x,
    cor = dd$estimate,
    p.value = dd$p.value
  )
}))

save(cor, file = "T_cor.Rds")
cor <- na.omit(cor)
cor <- cor[order(cor$cor, decreasing = T), ]
top <- 30
cor$label <- c(cor$gene[1:top], rep(NA, (nrow(cor) - top)))
# cor$p=ifelse(cor$cor>=0,2-(cor$p.value),(cor$p.value))
cor$p <- (nrow(cor):1)
p1 <- ggplot(data = cor) +
  geom_point(
    aes(
      x = cor,
      y = p
    ),
    alpha = 0.5,
    size = 1,
    color = "gray"
  ) +
  # geom_text_repel(aes(x=cor,y=p,label=label))+
  theme_classic() +
  geom_vline(xintercept = 0, lty = 5) +
  xlab("Correlation with granzyme expression") +
  ylab("Ranking") +
  annotate("rect",
    xmin = min(cor[1:top, "cor"]),
    xmax = Inf,
    ymin = 17000, ymax = Inf,
    alpha = 0.5,
    fill = "orange"
  )


data <- cor[1:top, ]
data$label <- factor(data$label, levels = rev(data$label))
p2 <- ggplot(data = data) +
  geom_segment(
    aes(
      xend = cor,
      x = 0,
      y = label,
      yend = label
    ),
    color = "grey"
  ) +
  geom_point(
    aes(
      x = cor,
      y = label,
      color = p.value
    ),
    size = 4
  ) +
  scale_color_gradient(high = "#FFE4B5", low = "orange") +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("R") +
  ylab(" ") +
  labs(title = paste("Top Gene: ", top))

p1 + p2 + plot_layout(widths = c(2, 1))

ggsave(filename = here("Fig/Fig5e3.pdf"), height = 6, width = 10)
