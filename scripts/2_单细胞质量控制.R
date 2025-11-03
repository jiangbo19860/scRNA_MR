rm(list = ls())

pacman::p_load(Seurat, dplyr, data.table, ggplot2, cowplot, patchwork, ggrepel, hdf5r, Matrix)

# 读取表达矩阵（注意使用 gzfile() 处理压缩文件）
counts <- readMM(gzfile("scRNA_read/matrix.mtx.gz"))
# 读取细胞条形码
barcodes <- readLines(gzfile("scRNA_read/barcodes.tsv.gz"))

# 读取基因信息
features <- read.table(gzfile("scRNA_read/features.tsv.gz"),
  header = FALSE,
  stringsAsFactors = FALSE
)

# 为基因表达矩阵添加行名和列名，使其成为一个“带注释的”表达矩阵。
rownames(counts) <- features[, 2] # 1表示基因ID（ENSG开头的）, 2表示基因名。
colnames(counts) <- barcodes # 细胞条形码

# 查看矩阵结构
dim(counts) # 输出：基因数 x 细胞数
counts[1:5, 1:5]
counts <- as(counts, "CsparseMatrix") # 要把counts转换为稀疏矩阵格式。
class(counts)

# 把稀疏矩阵转为Seurat对象 ---------------------------------------------------------
scRNA <- counts

# 去重，只保留第一次出现的基因
scRNA <- scRNA[!duplicated(rownames(scRNA)), ]

# 创建Seurat对象
scRNA <- CreateSeuratObject(
  counts = scRNA, # 输入的基因表达计数矩阵
  project = "scRNA", # 项目名称
  min.cells = 3, # 基因过滤阈值：仅保留在至少3个细胞中表达的基因
  min.features = 200 # 细胞过滤阈值：仅保留检测到至少200个基因的细胞
)

# 提取Seurat对象中的细胞元数据
meta <- scRNA@meta.data

colnames(scRNA@meta.data) # 打印元数据列名
str(scRNA@meta.data) # 单独打印元数据结构
head(scRNA@meta.data)

# 计算线粒体基因百分比（假设是人类数据，线粒体基因以"MT-"开头）
scRNA <- PercentageFeatureSet(scRNA, pattern = "^MT-", col.name = "percent.mt")

# 查看关键质量指标统计摘要
summary(scRNA@meta.data$nFeature_RNA) # 每个细胞检测到的基因数
summary(scRNA@meta.data$nCount_RNA) # 每个细胞的UMI总数
summary(scRNA@meta.data$percent.mt) # 线粒体基因比例

# 获取分组数量（用于配色）
group_num <- length(levels(scRNA@active.ident))
group_num

# 绘制QC指标小提琴图
# 在绘制小提琴图之前，添加以下代码计算中位数
median_feature <- median(scRNA@meta.data$nFeature_RNA, na.rm = TRUE) # 计算nFeature_RNA的中位数
median_count <- median(scRNA@meta.data$nCount_RNA, na.rm = TRUE) # 计算nCount_RNA的中位数
median_mt <- median(scRNA@meta.data$percent.mt, na.rm = TRUE) # 计算percent.mt的中位数

# 基本小提琴图（按活跃标识分组）
vln_plots <- VlnPlot(
  scRNA, # Seurat 对象，包含单细胞数据
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), # 指定要可视化的特征（基因数、UMI数、线粒体比例）
  cols = rainbow(group_num), # 使用彩虹色（rainbow）为每个分组生成颜色，group_num 是分组数量
  pt.size = 0.1, # 散点大小（每个点代表一个细胞，0.1 表示极小的点）
  ncol = 3 # 面板排列方式：一行显示3个图
)

vln_plots

# 为每个子图添加对应指标的中位数线
vln_plots[[1]] <- vln_plots[[1]] + # 第1个子图：nFeature_RNA
  geom_hline(yintercept = median_feature, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 1) # linewidth是线条粗细

vln_plots[[2]] <- vln_plots[[2]] + # 第2个子图：nCount_RNA
  geom_hline(yintercept = median_count, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 1)

vln_plots[[3]] <- vln_plots[[3]] + # 第3个子图：percent.mt
  geom_hline(yintercept = median_mt, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 1)

# 打印图形（此时无警告，中位数线正常显示）
print(vln_plots)

# 绘制基因数与UMI数散点图
feature_count_plot <- FeatureScatter(
  scRNA,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "orig.ident",
  cols = rainbow(group_num)
) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  ggtitle("nFeature_RNA vs nCount_RNA")

feature_count_plot


# 绘制基因数与线粒体比例散点图，正常应该为弱到中度的负相关，因为活性高的细胞转录活跃，表达的基因数多，nFeature_RNA高，且细胞膜完整，线粒体基因释放少，因此线粒体基因比例低。而活性低的细胞，转录活性下降，nFeature_RNA低，且细胞膜可能破裂，线粒体基因释放多，因此线粒体基因比例高。
feature_mt_plot <- FeatureScatter(
  scRNA,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt",
  group.by = "orig.ident",
  cols = rainbow(group_num)
) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  ggtitle("nFeatures vs percent.mt")
feature_mt_plot

# 组合散点图
combined_scatter <- plot_grid(feature_count_plot, feature_mt_plot, ncol = 2)
print("散点图组合:")
print(combined_scatter)

# 筛选高质量细胞,比较运算的结果是逻辑值（logical），
high_quality <- scRNA@meta.data$nFeature_RNA > 200 &
  scRNA@meta.data$nFeature_RNA < 5000 &
  scRNA@meta.data$percent.mt < 10

# 计算并打印高质量细胞比例，sum计算逻辑向量high_quality中 “TRUE 的数量”。length计算向量high_quality的长度，即包含的元素总个数。
hq_percentage <- sum(high_quality) / length(high_quality) * 100
round(hq_percentage, 2)

# 创建筛选后的Seurat对象
scRNA_hq <- subset(scRNA, cells = rownames(scRNA@meta.data[high_quality, ]))

# 定义带时间戳的日志输出函数（放在脚本开头）
log_info <- function(msg) {
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] INFO: ", msg, "\n", sep = "")
}

# 保存筛选后的高质量细胞对象
saveRDS(scRNA_hq, file = "scRNA_hq.rds")
log_info(paste("高质量细胞对象已保存至 scRNA_hq.rds，细胞数 =", ncol(scRNA_hq)))


# 可视化筛选前后的细胞分布
before_filter <- VlnPlot(
  scRNA,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  cols = "gray",
  pt.size = 0.1,
  ncol = 3
) + ggtitle("筛选前") +
  theme(
    axis.title.x = element_text(family = "ns_regular", size = 12), # 设置X轴标题字体
    axis.text.x = element_text(family = "ns_regular", size = 10) # 设置X轴刻度字体
  )

# 创建筛选后的临时对象
scRNA_hq <- subset(scRNA, cells = rownames(scRNA@meta.data[high_quality, ]))
after_filter <- VlnPlot(
  scRNA_hq,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  cols = "skyblue",
  pt.size = 0.1,
  ncol = 3
) + ggtitle("筛选后") +
  theme(
    axis.title.x = element_text(family = "ns_regular", size = 12), # 设置X轴标题字体
    axis.text.x = element_text(family = "ns_regular", size = 10) # 设置X轴刻度字体
  )

# 组合筛选前后的对比图
qc_comparison <- plot_grid(before_filter, after_filter, ncol = 1)
print("筛选前后对比:")
print(qc_comparison)
