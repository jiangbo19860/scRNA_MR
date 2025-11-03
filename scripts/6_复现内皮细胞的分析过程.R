## Fig2_1,592 endothelial cells, 6 clusters

# 清除当前工作空间所有对象，避免变量冲突
rm(list = ls())

# 加载所需R包（数据处理、可视化、单细胞分析工具等）
pacman::p_load(
  here, Seurat, SCopeLoomR, loomR, dplyr, ggplot2, cowplot, patchwork, ggrepel,
  SingleR, celldex, BiocParallel, scater, scran, SingleCellExperiment, data.table,
  hdf5r, Matrix, tidyr, stringr, forcats, pheatmap, RColorBrewer, viridis, ggridges, tibble
)

# 查看当前项目根目录路径
here()

# 0. 定义路径，加载数据 ------------------------------------------------------------
# 定义scRNA处理后数据的路径
rds_file <- here("1_data", "E-MTAB-6149", "scRNA_processed.Rds")
# 定义细胞视图数据的路径
ec_file <- here("1_data", "E-MTAB-6149", "EC.Cellview.Rds")

# 定义输出路径函数：所有结果保存到3_outputs/E-MTAB-6149_outputs目录
output_path <- function(...) {
  here("3_outputs", "E-MTAB-6149_outputs", ...)
}

# 若输出目录不存在则创建
if (!dir.exists(output_path())) {
  dir.create(output_path(), recursive = TRUE)
  message("已创建输出目录:", output_path())
}

# 读取细胞视图数据（含tsne坐标等）
load(ec_file)
str(tsne.data) # 查看t-SNE坐标数据的结构
unique(tsne.data$dbCluster) # 查看t-SNE数据中唯一的聚类标签
unique(tsne.data$V3) # 查看t-SNE数据中V3列的唯一值（此处均为NA）

# 读取scRNA处理后的Seurat对象
load(rds_file)
ls() # 查看当前工作空间的对象
scRNA # 查看Seurat对象的基本信息（基因数、细胞数等）
str(featuredata) # 查看基因特征注释数据的结构
str(tsne.data) # 再次确认t-SNE数据结构
str(scRNA@meta.data) # 查看Seurat对象的细胞元数据结构
scRNA@assays # 查看Seurat对象中的分析模块（如RNA）
names(scRNA@assays$RNA@layers) # 查看RNA模块中的数据层（如counts、data等）
str(scRNA@assays$RNA@layers$counts) # 查看原始表达计数矩阵的结构
dim(scRNA@assays$RNA@layers$counts) # 查看原始计数矩阵的维度（基因数×细胞数）

dim(scRNA) # 查看Seurat对象的维度（基因数×细胞数）= 22180×52698

class(log2cpm) # 查看log2cpm的类型（通常为矩阵，存储标准化后的表达量）

# 按log2cpm的列（细胞）筛选scRNA对象，只保留log2cpm列中的细胞，即分析内皮EC细胞
scRNA <- scRNA[, colnames(log2cpm)]

# 向scRNA添加t-SNE坐标数据作为元数据
scRNA <- AddMetaData(scRNA, metadata = tsne.data)
dim(scRNA) # 22180 × 1592，1592个细胞

# 对表达数据进行对数标准化（消除测序深度差异）
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# 筛选2000个高变基因（用于后续降维和聚类）
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10) # 取前10个高变基因

metadata <- scRNA@meta.data # 提取细胞元数据

## 数据中心化与回归
# 基于高变基因进行数据中心化（均值为0，方差为1）
scale.genes <- VariableFeatures(scRNA) # 用高变基因列表（减少内存占用）
scRNA <- ScaleData(scRNA, features = scale.genes)

# 基于高变基因进行PCA降维
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

# 绘制PCA图，按原始样本分组着色
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

# 绘制PCA肘图，确定合适的主成分数量，number of dimensions，要展示的主成分（PC）的数量（即维度数量）。
ElbowPlot(scRNA, ndims = 20, reduction = "pca")

pc.num <- 1:15 # 选择前15个主成分用于后续分析

# 基于选定的主成分计算细胞间距离
scRNA <- FindNeighbors(scRNA, dims = pc.num)

# 对细胞进行聚类（分辨率resolution=0.1，聚类数量较少,数字越大，分出的群越多）
scRNA <- FindClusters(scRNA, resolution = 0.1)

metadata <- scRNA@meta.data # 更新细胞元数据（包含聚类结果）

# 基于PCA结果进行t-SNE和UMAP降维（可视化用）
scRNA <- RunTSNE(scRNA, dims = pc.num)
scRNA <- RunUMAP(scRNA, dims = pc.num)

# 保存处理后的scRNA对象到输出目录
save(scRNA, file = output_path("scRNA_endo.Rds"))
# load(output_path("scRNA_endo.Rds")) # 重新加载保存的对象


# 1. 基于 t-SNE 降维的细胞聚类分布图_t-SNE visualization of endothelial cell clusters -----------------------------------------------------------
# 1.1 绘制t-SNE的基本图，按Seurat聚类结果着色，未指定 group.by 时，Seurat 默认使用Seurat聚类结果(seurat_clusters)列作为分组依据。
DimPlot(scRNA, reduction = "tsne")

# 绘制t-SNE 聚类可视化图，按dbCluster分组着色
DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne")

# 1.2 绘制基于dbCluster的t-SNE图（添加边框、隐藏图例、设置标题）并保存
DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.position = "none"
  ) +
  labs(title = "Endothelial cell")
# 保存图片到输出目录
ggsave(filename = output_path("Fig2a1_按细胞聚类绘制t-SNE图.pdf"), height = 8, width = 8)

# 1.3 绘制t-SNE图，按肿瘤来源（CellFromTumor）分组着色并美化
DimPlot(scRNA,
  reduction = "tsne",
  group.by = "CellFromTumor"
) +
  scale_color_manual(
    values = c("#2D5474", "#72C667"), # 先写哪个就在上面
    labels = c("Tumor", "Non-malignant")
  ) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.position = c(.01, .1)
  ) +
  labs(title = "Sample Origin")
# 保存图片到输出目录
ggsave(filename = output_path("Fig2a2_按细胞来源绘制t-SNE图.pdf"), height = 8, width = 8)

# 2. 绘制聚类中肿瘤/非肿瘤细胞的比例堆叠图 -----------------------------------------------------
# 准备聚类与样本来源的比例数据，赋值个tmp
tmp <- scRNA@meta.data # 1592行，18列

# 新增group列：区分肿瘤/非肿瘤细胞
tmp$group <- ifelse(tmp$CellFromTumor == "1", "Tumor", "Non-malignant")
# 调整聚类编号顺序（4到0）
tmp$seurat_clusters <- factor(tmp$seurat_clusters, levels = 4:0)

# 堆叠百分比柱状图（Stacked Percentage Bar Chart），展示细胞聚类（seurat_clusters）与样本分组（group）之间的比例关系。具体来说，它展示了每个细胞聚类中 “肿瘤细胞” 和 “非肿瘤细胞” 的比例分布。
p1 <- ggplot(
  tmp,
  aes(
    x = seurat_clusters,
    fill = group # 按照group列（肿瘤/非肿瘤）进行填充
  )
) +
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  scale_fill_manual(values = c("#72C667", "#2D5474")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(" ") +
  ylab(" ") +
  theme(legend.position = "top") +
  theme(
    axis.line.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), axis.text.y = element_text(size = 14)
  ) +
  coord_flip()
p1
# 绘制聚类中不同患者的比例堆叠图
p2 <- ggplot(
  tmp,
  aes(x = seurat_clusters, fill = PatientNumber)
) +
  geom_bar(stat = "count", position = "fill") +
  theme_classic() +
  scale_fill_manual(values = c("#268A24", "#8EEE8B", "#FB6346", "#FBD51A", "#28507D")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(" ") +
  ylab(" ") +
  theme(legend.position = "top") +
  theme(
    axis.line.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), axis.text.y = element_text(size = 14)
  ) +
  coord_flip()
p2

# 组合两个图并保存
p1 + p2 + plot_layout(ncol = 2)
ggsave(output_path("Fig1c_endo_堆叠图.pdf"), height = 4, width = 10)


################################################################################
# 绘制指定基因（FLT1/PDPN等）在t-SNE上的表达分布图
p <- FeaturePlot(scRNA,
  reduction = "tsne",
  features = c("FLT1", "PDPN", "MT2A", "HSPG2"),
  ncol = 2, # 2列排列
  cols = c("gray", "red") # 表达量从低到高：灰色到红色
)
# 美化图（添加边框、隐藏图例）
p & theme(
  panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
  legend.position = "none"
)
# 保存图片到输出目录
ggsave(output_path("Fig2b.pdf"), height = 8, width = 8)

######################################################################################
# 加载通路分析相关包
library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

# 读取Hallmark基因集（.gmt格式）
gs <- read.gmt("h.all.v7.4.symbols.gmt")
# 去除基因集名称中的"HALLMARK_"前缀
gs$term <- gsub("HALLMARK_", "", gs$term)

# 转换为基因集列表（用于GSVA分析）
gs.list <- gs %>%
  split(.$term) %>%
  lapply("[[", 2)

# 进行GSVA分析（计算每个细胞的通路富集分数）
gsva_es <- gsva(as.matrix(scRNA@assays$RNA@counts),
  gs.list,
  method = "ssgsea", # 采用ssGSEA算法
  abs.ranking = T
)

# 为scRNA对象添加肿瘤/非肿瘤分组信息
scRNA$group <- ifelse(scRNA$CellFromTumor == "1", "tumor", "Non_malignant")

# 构建分组矩阵（用于limma差异分析）
group_list <- data.frame(
  sample = colnames(gsva_es),
  group = scRNA$group
)
head(group_list)

# 构建设计矩阵
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

# 构建对比矩阵（肿瘤vs非肿瘤）
contrast.matrix <- makeContrasts(tumor - Non_malignant, levels = design)

# 用limma进行通路差异分析
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取所有通路的差异结果（按t值排序）
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

# 保存差异通路结果到输出目录
write.csv(x, output_path("endo_gsva_limma.csv"), quote = F)

# 整理通路t值数据用于绘图
df <- data.frame(ID = rownames(x), score = x$t)
# 按t值分组（高/中/低）
cutoff <- 2
df$group <- cut(df$score,
  breaks = c(-Inf, -cutoff, cutoff, Inf),
  labels = c(1, 2, 3)
)
# 按t值排序
sortdf <- df[order(df$score), ]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

# 绘制通路t值条形图（按t值排序）
ggplot(sortdf, aes(ID, score, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("palegreen3", "snow3", "dodgerblue4"), guide = FALSE) +
  geom_hline(yintercept = c(-cutoff, cutoff), color = "white", linetype = 2, size = 0.3) +
  geom_text(
    data = subset(df, score < 0), aes(x = ID, y = 0.1, label = ID, color = group),
    size = 3, hjust = 0
  ) +
  geom_text(
    data = subset(df, score > 0), aes(x = ID, y = -0.1, label = ID, color = group),
    size = 3, hjust = 1
  ) +
  scale_colour_manual(values = c("black", "snow3", "black"), guide = FALSE) +
  xlab("") +
  ylab("t value of GSVA score, tumor versus non-malignant") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 0.6)) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
# 保存图片到输出目录
ggsave(filename = output_path("Fig2c.pdf"), height = 10, width = 7)

#################################################################################
# 定义关注的基因列表
gene_name <- c("HSPG2", "ANGPT2", "HIF1A", "MMP2", "CTGF", "NOTCH1")

# 提取基因表达量并关联聚类信息
data_tmp <- data.frame(
  group = scRNA$seurat_clusters,
  row.names = rownames(scRNA@meta.data)
)
# 循环提取每个基因的表达量
for (i in 1:6) {
  tmp <- FeaturePlot(scRNA, features = gene_name[i]) # 获取基因表达数据
  tmp <- tmp$data
  tmp <- data.frame(gene = as.numeric(tmp[, 4]), row.names = rownames(scRNA@meta.data))
  colnames(tmp) <- gene_name[i]
  data_tmp <- cbind(tmp, data_tmp) # 合并到数据框
}

# 筛选特定聚类（0-3）的数据
data_tmp <- subset(data_tmp, group %in% c(0, 1, 2, 3))
data_tmp$group <- factor(data_tmp$group, levels = c(0, 1, 2, 3)) # 调整聚类顺序

# 转换数据格式用于小提琴图
new_vio <- reshape2::melt(data = data_tmp, value.name = "exp")
new_vio$fill <- ifelse(new_vio$group %in% c(0, 1), "A", "B") # 分组填充颜色

# 绘制小提琴图（按基因分面）
ggplot(new_vio, aes(x = group, y = exp)) +
  geom_violin(aes(fill = fill), show.legend = F) +
  scale_fill_manual(values = c("#7CCD7C", "#36648B")) +
  theme_bw() +
  facet_grid(variable ~ .) + # 按基因纵向分面
  xlab("") +
  ylab("") +
  theme(
    panel.grid = element_blank(), strip.background.x = element_blank(),
    panel.border = element_rect(size = 1.2), axis.line = element_line(size = 1.2),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 10)
  )

#################################################################################
# 加载堆叠小提琴图函数
source("stackvlion.R")
# 绘制堆叠小提琴图（展示基因在不同聚类的表达）
StackedVlnPlot(scRNA,
  features = gene_name, idents = c(0, 1, 2, 3, 4),
  cols = c("#7CCD7C", "#7CCD7C", "#36648B", "#36648B", "#36648B")
)
# 保存图片到输出目录
ggsave(filename = output_path("Fig2f.pdf"), height = 5, width = 4)
