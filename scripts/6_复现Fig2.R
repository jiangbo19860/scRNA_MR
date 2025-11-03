rm(list = ls())
pacman::p_load(
  here, Seurat, SCopeLoomR, loomR, dplyr, ggplot2, cowplot, patchwork, ggrepel, SingleR, celldex, BiocParallel, scater, scran, SingleCellExperiment, data.table, hdf5r, Matrix, tidyr, stringr, forcats, pheatmap, RColorBrewer, viridis, ggridges, tibble
)
here()

# 0. 定义路径，加载数据 ------------------------------------------------------------
rds_file <- here("1_data", "E-MTAB-6149", "scRNA_processed.Rds")
ec_file <- here("1_data", "E-MTAB-6149", "EC.Cellview.Rds")

# 定义输出路径函数，所有输出文件保存在3_outputs目录下的E-MTAB-6149_outputs子目录中。
output_path <- function(...) {
  here("3_outputs", "E-MTAB-6149_outputs", ...)
}

# 创建输出目录（如果不存在）
if (!dir.exists(output_path())) {
  dir.create(output_path(), recursive = TRUE)
  message("已创建输出目录:", output_path())
}

# 读取Rds文件
load(ec_file) # 读取Rds文件，加载数据到当前工作空间
str(tsne.data) # 查看tsne.data的结构
unique(tsne.data$dbCluster) # 查看tsne.data中细胞的唯一标识符
unique(tsne.data$V3) # 查看tsne.data中细胞的唯一标识符)
load(rds_file) # 读取Rds文件，加载数据到当前工作空间
ls()
scRNA
str(featuredata)
str(tsne.data)
str(scRNA@meta.data)
scRNA@assays # 查看Seurat对象中的assays（数据层）
names(scRNA@assays$RNA@layers) # 查看RNA模块中所有层的结构
str(scRNA@assays$RNA@layers$counts) # 查看RNA模块中counts层（原始表达计数）的结构
dim(scRNA@assays$RNA@layers$counts) # 查看RNA模块中counts层的维度


dim(scRNA) # 查看scRNA对象的维度

class(log2cpm) # log2-transformed Counts Per Million”（对数转换后的每百万次测序读数）。

# 按列（细胞）筛选scRNA这个Seurat对象，只保留列名在colnames(log2cpm)中的列（细胞条形码-唯一标识符）。
scRNA <- scRNA[, colnames(log2cpm)]

# 向scRNA 这个 Seurat 对象中添加新的元数据（细胞水平的信息），这里添加的是 tsne.data（通常是提前计算好的 t-SNE 降维坐标）。
scRNA <- AddMetaData(scRNA, metadata = tsne.data)

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)

metadata <- scRNA@meta.data
## 数据中心化与回归
# 数据中心化
# scale.genes <-  rownames(scRNA)
scale.genes <- VariableFeatures(scRNA) # 内存不够可用此命令代替上一行
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

save(scRNA, file = "scRNA_endo.Rds")
load("scRNA_endo.Rds")

DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne")
DimPlot(scRNA, reduction = "tsne")

DimPlot(scRNA, reduction = "tsne") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    legend.position = "none"
  ) +
  labs(title = "Endothelial cell")

ggsave(filename = "Fig/Fig2a.pdf", height = 8, width = 8)

####################################################################################

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

ggsave(filename = "Fig/Fig2a1.pdf", height = 8, width = 8)


####################################################################################
tmp <- scRNA@meta.data
tmp$group <- ifelse(tmp$CellFromTumor == "1", "Tumor", "Non-malignant")
tmp$seurat_clusters <- factor(tmp$seurat_clusters, levels = 4:0)

p1 <- ggplot(
  tmp,
  aes(
    x = seurat_clusters,
    fill = group
  )
) +
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  scale_fill_manual(values = c("#72C667", "#2D5474")) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5
  )) +
  xlab(" ") +
  ylab(" ") +
  theme(legend.position = "top") +
  theme(axis.line.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  coord_flip()

p2 <- ggplot(
  tmp,
  aes(x = seurat_clusters, fill = PatientNumber)
) +
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  scale_fill_manual(values = c(
    "#268A24",
    "#8EEE8B",
    "#FB6346",
    "#FBD51A",
    "#28507D"
  )) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5
  )) +
  xlab(" ") +
  ylab(" ") +
  theme(legend.position = "top") +
  theme(axis.line.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  coord_flip()
p1 + p2 + plot_layout(ncol = 2)

ggsave("Fig/Fig2_endo_bar.pdf", height = 4, width = 10)


################################################################################

p <- FeaturePlot(scRNA,
  reduction = "tsne",
  features = c(
    "FLT1",
    "PDPN",
    "MT2A",
    "HSPG2"
  ),
  ncol = 2,
  cols = c("gray", "red")
)
p & theme(
  panel.border = element_rect(
    fill = NA, color = "black",
    size = 1, linetype = "solid"
  ),
  legend.position = "none"
)

ggsave("Fig/Fig2b.pdf", height = 8, width = 8)

######################################################################################

library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

gs <- read.gmt("h.all.v7.4.symbols.gmt")

gs$term <- gsub("HALLMARK_", "", gs$term)

gs.list <- gs %>%
  split(.$term) %>%
  lapply("[[", 2)

gsva_es <- gsva(as.matrix(scRNA@assays$RNA@counts),
  gs.list,
  method = "ssgsea",
  abs.ranking = T
)

scRNA$group <- ifelse(scRNA$CellFromTumor == "1", "tumor", "Non_malignant")

group_list <- data.frame(
  sample = colnames(gsva_es),
  group = scRNA$group
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

x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

head(x)

write.csv(x, "endo_gsva_limma.csv", quote = F)

df <- data.frame(ID = rownames(x), score = x$t)

# 按照score的值分组
cutoff <- 2
df$group <- cut(df$score,
  breaks = c(-Inf, -cutoff, cutoff, Inf),
  labels = c(1, 2, 3)
)

# 按照score排序
sortdf <- df[order(df$score), ]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID,
  score,
  fill = group
)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "palegreen3",
      "snow3",
      "dodgerblue4"
    ),
    guide = FALSE
  ) +
  geom_hline(
    yintercept = c(-cutoff, cutoff),
    color = "white",
    linetype = 2, # 画虚线
    size = 0.3
  ) +
  geom_text(
    data = subset(df, score < 0),
    aes(
      x = ID,
      y = 0.1,
      label = ID,
      color = group
    ), # bar跟坐标轴间留出间隙
    size = 3, # 字的大小
    hjust = 0
  ) + # 字的对齐方式
  geom_text(
    data = subset(df, score > 0),
    aes(
      x = ID,
      y = -0.1,
      label = ID,
      color = group
    ),
    size = 3,
    hjust = 1
  ) +
  scale_colour_manual(values = c("black", "snow3", "black"), guide = FALSE) +
  xlab("") +
  ylab("t value of GSVA score, tumor \n versus non-malignant") +
  theme_bw() + # 去除背景色
  theme(panel.grid = element_blank()) + # 去除网格线
  theme(panel.border = element_rect(size = 0.6)) + # 边框粗细
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) # 去除y轴

ggsave(filename = "Fig/Fig2c.pdf", height = 10, width = 7)

#################################################################################

gene_name <- c(
  "HSPG2",
  "ANGPT2",
  "HIF1A",
  "MMP2",
  "CTGF",
  "NOTCH1"
)

data_tmp <- data.frame(
  group = scRNA$seurat_clusters,
  row.names = rownames(scRNA@meta.data)
)

for (i in 1:6) {
  tmp <- FeaturePlot(scRNA, features = gene_name[i])
  tmp <- tmp$data
  tmp <- data.frame(
    gene = as.numeric(tmp[, 4]),
    row.names = rownames(scRNA@meta.data)
  )
  colnames(tmp) <- gene_name[i]
  data_tmp <- cbind(tmp, data_tmp)
}

data_tmp <- subset(data_tmp, group %in% c(0, 1, 2, 3))

data_tmp$group <- factor(data_tmp$group, levels = c(0, 1, 2, 3))

new_vio <- reshape2::melt(
  data = data_tmp,
  value.name = "exp"
)
new_vio$fill <- ifelse(new_vio$group %in% c(0, 1), "A", "B")

ggplot(new_vio, aes(
  x = group,
  y = exp
)) +
  geom_violin(aes(fill = fill),
    show.legend = F
  ) +
  scale_fill_manual(values = c("#7CCD7C", "#36648B")) +
  theme_bw() +
  facet_grid(variable ~ .) +
  xlab("") +
  ylab("") +
  theme(
    panel.grid = element_blank(),
    strip.background.x = element_blank(),
    panel.border = element_rect(size = 1.2),
    axis.line = element_line(size = 1.2),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 10)
  )

#################################################################################

source("stackvlion.R")
StackedVlnPlot(scRNA,
  features = gene_name, idents = c(0, 1, 2, 3, 4),
  cols = c(
    "#7CCD7C", "#7CCD7C",
    "#36648B", "#36648B", "#36648B"
  )
)

ggsave(filename = "Fig/Fig2f.pdf", height = 5, width = 4)
