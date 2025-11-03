### 复现Fig3_纤维母细胞的结果

rm(list = ls())
pacman::p_load(
  Seurat, # 单细胞数据分析核心包
  dplyr, # 数据处理
  limma, # 差异分析
  stringr, # 字符串处理
  GSVA, # 基因集变异分析
  clusterProfiler, # 富集分析
  data.table, # 高效数据框操作
  ggplot2, # 可视化
  cowplot, patchwork, # 图形组合
  ggrepel, # 标签避免重叠
  hdf5r, # 读取HDF5格式文件
  BiocParallel # 并行计算
)
here() # 确认当前项目路径


# 数据加载，筛选细胞，添加元数据，标准化，高变基因选择，归一化ScaleData,PCA,ElbowPlot, FindNeighbors, FindClusters, RunTSNE/RunUMAP,  保存为.Rds--------------------------------------------------------------------
# 加载纤维母细胞相关数据和处理后的单细胞数据
load(here("1_data", "E-MTAB-6149", "Fibro.Cellview.Rds"))
load(here("1_data", "E-MTAB-6149", "scRNA_processed.Rds"))
# 查看log2CPM标准化后的表达矩阵（前3行3列）
log2cpm[1:3, 1:3]
# 提取Seurat对象的元数据, 每行代表一个细胞。
meta <- scRNA@meta.data
counts <- scRNA@assays$RNA@layers$counts # 原始计数矩阵（基因×细胞）
data <- scRNA@assays$RNA@layers$data # 标准化后的表达矩阵（基因×细胞）
dim(counts) # 22180 genes x 52698 cells
dim(data) # 22180 genes x 52698 cells
counts[1:5, 1:5]

# 根据log2cpm矩阵的列名（细胞ID）筛选Seurat对象中的细胞（保留匹配的细胞）
scRNA <- scRNA[, colnames(log2cpm)]
# 向Seurat对象添加t-SNE坐标元数据（用于后续可视化）
scRNA <- AddMetaData(scRNA, metadata = tsne.data)
# 对表达数据进行标准化（LogNormalize方法：log2(计数/总和 + 1)，缩放因子10000）
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
# 选择高变基因（使用vst方法，筛选2000个变异最大的基因）
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# 提取排名前10的高变基因（查看筛选结果）
top10 <- head(VariableFeatures(scRNA), 10)

# 对表达矩阵进行中心化（均值为0，标准差为1）
scRNA <- ScaleData(scRNA)
# 基于高变基因进行主成分分析（PCA）
scRNA <- RunPCA(scRNA)
# 可视化PCA结果（按原始样本来源分组，查看是否存在批次效应）
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
# 绘制Elbow图（碎石图），确定用于后续分析的最佳主成分（PC）数量
# 拐点处的PC数量即包含大部分变异信息（此处根据图选择前15个PC）
ElbowPlot(scRNA, ndims = 20, reduction = "pca")
# 设置要使用的主成分数量（1到15）
pc.num <- 1:15
# 基于PCA结果构建细胞近邻图（用于聚类）
scRNA <- FindNeighbors(scRNA, dims = pc.num) # 使用前15个PC
# 对细胞进行聚类（resolution=0.2：聚类分辨率，值越小簇数量越少）
scRNA <- FindClusters(scRNA, resolution = 0.2) # 聚类结果保存在seurat_clusters列中
# 提取聚类后的元数据（包含每个细胞的簇标签）
metadata <- scRNA@meta.data
# 进行t-SNE降维（将高维PCA结果转换为2D可视化）
scRNA <- RunTSNE(scRNA, dims = pc.num)
# 进行UMAP降维（另一种常用的2D可视化方法，保留局部结构更好）
scRNA <- RunUMAP(scRNA, dims = pc.num)
# 保存处理后的Seurat对象（包含t-SNE和UMAP结果），方便后续分析，避免重复计算
save(scRNA, file = "scRNA_fibro_0718.Rds")

unique(scRNA$seurat_clusters) # 查看聚类结果（簇标签）
unique(scRNA$CellFromTumor) # 查看细胞来源（肿瘤或非肿瘤）
unique(scRNA$dbCluster) # 查看自定义簇标签（dbCluster）1-7
unique(scRNA$PatientNumber) # 查看患者编号（1到5）

# Fig1a_按dbCluster分组绘制t-SNE图（查看自定义簇标签分布）------------------------------------------------------------------
# load("scRNA_fibro_0718.Rds")

DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne")
# 绘制t-SNE图（默认按FindClusters的结果分组，展示聚类效果）
DimPlot(scRNA, reduction = "tsne")

# 绘制按FindClusters函数生成的seurat_clusters分组的t-SNE图
DimPlot(scRNA, group.by = "seurat_clusters", reduction = "tsne", label = TRUE) +
  theme(
    panel.border = element_rect(
      fill = NA, color = "black",
      size = 1, linetype = "solid"
    ),
    legend.position = "none"
  ) +
  labs(title = "Fibroblasts")

ggsave(filename = "Fig/Fig3a1_按FindClusters函数生成的seurat_clusters分组的t-SNE图.pdf", height = 8, width = 8)

# 绘制按dbCluster分组的t-SNE图
DimPlot(scRNA, group.by = "dbCluster", reduction = "tsne", label = TRUE) +
  theme(
    panel.border = element_rect(
      fill = NA, color = "black",
      size = 1, linetype = "solid"
    ),
    legend.position = "none"
  ) +
  labs(title = "Fibroblasts")

ggsave(filename = "Fig/Fig3a1_按dbCluster分组的t-SNE图.pdf", height = 8, width = 8)

####################################################################################

# Fig3a2_DimPlot展示在t-SNE降维空间上细胞的来源分组（肿瘤或非肿瘤） -------------------------------

# 在 t-SNE 降维空间上按细胞来源（肿瘤/非肿瘤）着色并可视化
DimPlot(
  scRNA, # Seurat 对象，包含单细胞数据
  reduction = "tsne", # 使用 t-SNE 降维结果
  group.by = "CellFromTumor" # 按 CellFromTumor 列分组（值为 1 或 0）
) +
  # 手动设置颜色映射：蓝色(#2D5474)表示肿瘤，绿色(#72C667)表示非肿瘤
  scale_color_manual(
    values = c("#2D5474", "#72C667"),
    labels = c("Tumor", "Non-malignant") # 图例标签
  ) +
  # 自定义图形主题
  theme(
    panel.border = element_rect( # 添加黑色边框
      fill = NA, color = "black",
      size = 1, linetype = "solid"
    ),
    legend.position = c(.01, .1) # 图例位置：左下角(0.01, 0.1)
  ) +
  # 设置图形标题
  labs(title = "Sample Origin")

# 保存图形为 PDF 文件（8×8 英寸）
ggsave(filename = "Fig/Fig3a2_DimPlot展示在t-SNE降维空间上细胞来源于肿瘤还是非肿瘤.pdf", height = 8, width = 8)


####################################################################################

# Fig1d_4_ggplot2包绘制了两个水平堆叠柱状图，分别展示细胞簇的肿瘤/非肿瘤来源比例和患者来源比例--------
# 提取Seurat对象的元数据
tmp <- scRNA@meta.data

# 创建新列"group"，根据CellFromTumor的值标记为"Tumor"或"Non-malignant"
tmp$group <- ifelse(tmp$CellFromTumor == "1", "Tumor", "Non-malignant")

# 将seurat_clusters转换为因子，并指定顺序为5到0（从顶部到底部排列）
tmp$seurat_clusters <- factor(tmp$seurat_clusters, levels = 5:0)

# 细胞簇的肿瘤 / 非肿瘤来源比例
p1 <- ggplot(
  tmp, # 数据源
  aes(
    x = seurat_clusters, # x轴为细胞簇（实际水平显示，因coord_flip）
    fill = group # 填充色为肿瘤/非肿瘤来源
  )
) +
  # 统计每个簇中各组的细胞数，并按比例填充（position="fill"）
  geom_bar(stat = "count", position = "fill") +
  # 使用经典主题（无网格背景）
  theme_classic() +
  # 手动指定填充色：绿色=非肿瘤，蓝色=肿瘤
  scale_fill_manual(values = c("#72C667", "#2D5474")) +
  # 隐藏/调整x轴元素（因水平翻转，实际为垂直轴）
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(" ") + # 移除x轴标题
  theme(axis.line.x = element_blank()) + # 移除x轴线
  theme(axis.text.x = element_blank()) + # 移除x轴文本
  theme(axis.ticks.x = element_blank()) + # 移除x轴刻度

  # 隐藏/调整y轴元素（因水平翻转，实际为水平轴）
  ylab(" ") + # 移除y轴标题
  theme(axis.text.y = element_text(size = 14)) + # 设置y轴文本大小
  theme(axis.line.y = element_blank()) + # 移除y轴线
  theme(axis.ticks.y = element_blank()) + # 移除y轴刻度

  # 将图例放在顶部
  theme(legend.position = "top") +
  # 水平翻转坐标轴（使柱状图水平显示）
  coord_flip()

# 细胞簇的患者来源比例
p2 <- ggplot(
  tmp, # 数据源
  aes(
    x = seurat_clusters, # x轴为细胞簇（实际水平显示，因coord_flip）
    fill = PatientNumber # 填充色为患者编号
  )
) +
  # 统计每个簇中各组的细胞数，并按比例填充（position="fill"）
  geom_bar(stat = "count", position = "fill") +
  # 使用经典主题（无网格背景）
  theme_classic() +
  # 手动指定填充色：不同患者使用不同颜色
  scale_fill_manual(values = c(
    "#268A24", "#8EEE8B", "#FB6346", "#FBD51A", "#28507D"
  )) +
  # 隐藏/调整x轴元素（因水平翻转，实际为垂直轴）
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(" ") + # 移除x轴标题
  theme(axis.line.x = element_blank()) + # 移除x轴线
  theme(axis.text.x = element_blank()) + # 移除x轴文本
  theme(axis.ticks.x = element_blank()) + # 移除x轴刻度

  # 隐藏/调整y轴元素（因水平翻转，实际为水平轴）
  ylab(" ") + # 移除y轴标题
  theme(axis.text.y = element_text(size = 14)) + # 设置y轴文本大小
  theme(axis.line.y = element_blank()) + # 移除y轴线
  theme(axis.ticks.y = element_blank()) + # 移除y轴刻度

  # 将图例放在顶部
  theme(legend.position = "top") +
  # 水平翻转坐标轴（使柱状图水平显示）
  coord_flip()

# 使用patchwork包将两个图形并排显示
p_combined <- p1 + p2 + plot_layout(ncol = 2)

# 正确图片
ggsave(
  filename = "Fig/Fig1d_4_水平堆叠柱状图展示纤维母细胞簇的肿瘤和非肿瘤来源比例和患者来源比例.pdf",
  plot = p_combined, # 你的组合图对象（如p1 + p2）
  height = 4,
  width = 10
)


# Fig3b_用FeaturePlot(散点图)在t-SNE降维空间上可视化多个基因的表达分布 -------------------------------------------
# 绘制基因表达的特征图（FeaturePlot）
p <- FeaturePlot(
  scRNA, # Seurat对象，包含单细胞数据
  reduction = "tsne", # 使用t-SNE降维结果作为可视化空间
  features = c( # 指定要可视化的基因列表
    "COL10A1", "COL4A1", # 胶原蛋白相关基因
    "PLA2G2A", "MMP3", # 酶类相关基因
    "FIGF", "CCL2" # 生长因子和趋化因子
  ),
  ncol = 3, # 每行显示3个基因的表达图
  cols = c("gray", "red") # 颜色映射：低表达为灰色，高表达为红色
)

# 修改图形主题：添加黑色边框，移除图例
p & theme(
  panel.border = element_rect( # 设置每个子图的边框
    fill = NA, # 边框内部不填充颜色
    color = "black", # 边框颜色为黑色
    size = 1, # 边框粗细为1磅
    linetype = "solid" # 边框线型为实线
  ),
  legend.position = "none" # 移除图例（颜色条）
)

# 保存图形为PDF文件
ggsave("Fig/Fig3b_FeaturePlot散点图展示多个基因在t-SNE降维空间上的表达分布.pdf", height = 6, width = 9)

# Fig3e_热图_基因集富集分析（GSEA）后不同细胞簇（clusters）中各通路的活性差异。 ----------------------------
# 基因集富集分析（GSEA），热图可视化不同细胞簇（clusters）中各通路的活性差异。
library(GSVA) # 执行基因集变异分析（GSVA）和ssGSEA
library(clusterProfiler) # 基因富集分析和可视化
library(limma) # 差异表达分析和线性模型
library(stringr) # 字符串处理
library(ggplot2) # 高级绘图

# 读取Hallmark基因集（GMT格式）
gs <- read.gmt(here("1_data", "GSEA", "h.all.v7.4.symbols.gmt"))

# 去除通路名称前缀"HALLMARK_"（如将"HALLMARK_APOPTOSIS"改为"APOPTOSIS"）
gs$term <- gsub("HALLMARK_", "", gs$term)

# 将数据框转换为通路-基因列表格式（每个通路对应一个基因向量）。gs.list是一个列表，每个元素对应一个通路及其包含的基因。
gs.list <- gs %>%
  split(.$term) %>% # 按通路名称分组
  lapply("[[", 2) # 提取每个分组中的第二列（基因名称列）

# 从Seurat对象中提取原始计数矩阵（Seurat 5+兼容）
counts_matrix <- GetAssayData(
  object = scRNA,
  assay = "RNA", # 指定assay（通常为RNA）
  layer = "counts" # 指定层为原始计数（对应旧版的@counts）
)

counts_matrix <- as.matrix(counts_matrix)

# 创建ssGSEA参数对象
ssgsea_param <- ssgseaParam(
  expr = counts_matrix, # 输入：标准化后的表达矩阵（基因×细胞）
  geneSets = gs.list # 输入：基因集列表（50个核心通路）
)
# 设置并行计算（3个CPU核心）
register(MulticoreParam(workers = 3)) # 3个核心，可根据电脑配置调整

# 重要：运行ssGSEA分析, 使用参数对象调用gsva(), 得到的gsva_result是一个通路 × 细胞的矩阵（基因集富集分数矩阵)。行：信号通路（如ADIPOGENESIS、APOPTOSIS等）；列: 细胞（与输入表达矩阵的细胞一一对应）；值：通路富集分数（数值越大，该通路在对应细胞中的活性越高）。
gsva_result <- gsva(
  param = ssgsea_param, # 传入参数对象（包含数据和算法设置）
  BPPARAM = bpparam(), # 应用并行计算参数（使用3个核心）
  verbose = FALSE # 不输出中间过程（静默运行）
)
class(gsva_result)

rownames(gsva_result)
head(rownames(scRNA@meta.data))

## 按细胞簇（clusters）汇总通路活性

# 获取scRNA中实际存在的簇编号（自动识别，无需硬编码）
existing_clusters <- sort(unique(scRNA$seurat_clusters))
cat("实际存在的簇：", existing_clusters, "\n")


# 初始化一个空的数据框 data_gs（行=通路，列=簇），专门用于存储后续分析结果。创建一个没有任何列的数据框（data.frame()）。将行名设置为 gsva_result 的行名（即 50 个 Hallmark 通路的名称，如 "APOPTOSIS"、"DNA_REPAIR" 等）。
data_gs <- data.frame(row.names = rownames(gsva_result)) # 初始化为空数据框（仅行名）

# 遍历实际存在的簇（existing_clusters是通过unique(scRNA$seurat_clusters)获取的非重复簇编号向量）
for (x in existing_clusters) {
  # 提取当前簇x的所有细胞ID：
  # - scRNA@meta.data存储了细胞的元数据（如聚类结果），每一行对应一个细胞
  # - which(scRNA$seurat_clusters == x)返回属于簇x的细胞在meta.data中的索引
  # - rownames(...)将索引转换为对应的细胞ID（如"AACGGTACCTTCGC_1"）
  group <- rownames(scRNA@meta.data)[which(scRNA$seurat_clusters == x)]

  # 检查该簇是否为空（即没有细胞被分配到该簇）：
  # - length(group) == 0表示该簇没有细胞
  # - continue语句跳过当前循环，继续处理下一个簇
  if (length(group) == 0) {
    cat("警告：簇", x, "为空，跳过\n")
    continue
  }

  # 提取属于当前簇的所有细胞的通路活性数据：
  # - gsva_result是通路×细胞的矩阵（50行×细胞数列）
  # - gsva_result[, group]表示提取所有行（通路）中，列名为group（即簇x的细胞ID）的子集
  gsva_es_tmp <- gsva_result[, group]

  # 计算每个通路在当前簇中的平均活性：
  # - apply(..., 1, mean)对矩阵的每一行（通路）计算均值
  # - 结果是一个长度为50的向量，每个元素对应一个通路的平均活性
  gsva_es_tmp <- apply(gsva_es_tmp, 1, mean)

  # 将向量转换为数据框格式：
  # - 便于后续添加列名和合并操作
  # - 转换后的数据框有1列，行数为50（对应50个通路）
  gsva_es_tmp <- as.data.frame(gsva_es_tmp)

  # 设置数据框的列名为当前簇编号：
  # - 例如，如果x=1，则列名设为"1"
  # - 这一步确保最终结果中每列的名称对应正确的簇
  colnames(gsva_es_tmp) <- x

  # 将当前簇的通路活性均值添加到结果数据框：
  # - cbind()函数按列合并数据框
  # - 每次循环添加一列，最终data_gs的列数等于实际存在的簇数
  data_gs <- cbind(data_gs, gsva_es_tmp)
}

# 确认结果数据框的列顺序与existing_clusters一致
cat("数据框列名顺序：", colnames(data_gs), "\n")

dim(data_gs) # 应输出：50 6（50个通路 × 6个cluster）

# 绘制热图（使用完整的data_gs）
pdf("Fig/Fig3e_heatmap展示gsva分析后不同通路在不同簇细胞中的活性score.pdf", height = 10, width = 8)
pheatmap(
  data_gs,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_cols = FALSE, # 不聚类列
  cluster_rows = FALSE, # 不聚类行
  scale = "row", # 对数据按行（通路）进行标准化处理。
  display_numbers = FALSE, # 不在热图的每个单元格中显示具体数值。
  fontsize_col = 8, # 设置热图 x 轴（列名，即簇编号）的字体大小为 8 号。
  main = "Pathway Activity Heatmap" # 添加标题
)
dev.off()



# 另一种方法绘制Fig3e_热图展示GSEA的结果 ------------------------------------------------------
# 获取实际存在的簇（替换0:5）
existing_clusters <- sort(unique(scRNA$seurat_clusters))
# 初始化data（行=通路）
data <- data.frame(row.names = rownames(gsva_result))

# 只循环实际存在的簇
for (i in existing_clusters) {
  # 标记目标簇为"A"，其他为"B"
  scRNA$group <- ifelse(scRNA$seurat_clusters == i, "A", "B")

  # 确保样本名与gsva_result的列名一致（细胞ID顺序必须匹配）
  group_list <- data.frame(
    sample = colnames(gsva_result), # gsva_result的列是细胞ID
    group = scRNA$group[match(colnames(gsva_result), rownames(scRNA@meta.data))] # 按细胞ID匹配group
  )

  # 检查A组是否有细胞（避免空组）
  if (sum(group_list$group == "A") == 0) {
    cat("警告：簇", i, "无细胞，跳过\n")
    next
  }

  # 差异分析（limma）
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- group_list$sample
  contrast.matrix <- makeContrasts(A - B, levels = design)
  fit <- lmFit(gsva_result, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  # 提取所有通路的t值（n=Inf确保不过滤）
  x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

  # 确保x包含所有通路（补充缺失通路的t值为NA）
  all_paths <- rownames(gsva_result)
  x_matched <- data.frame(
    t = x[match(all_paths, rownames(x)), "t"], # 按通路名匹配t值
    row.names = all_paths
  )

  # 合并到data
  colnames(x_matched) <- as.character(i)
  data <- cbind(data, x_matched)
}
# 检查是否有有效列
dim(data) # 应输出：通路数 × 实际簇数（如50 × 6）
# 检查是否有NA（可选）
sum(is.na(data)) # 若过多，需排查差异分析问题
class(data) # 确保是数据框

library(pheatmap)

pdf("Fig/Fig3e_gsva结果heatmap的另一种绘制方法.pdf", width = 10, height = 12)
pheatmap(
  data,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50), # 颜色方案（从蓝到白到红）
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "none" # 不标准化，直接用原始值
)
dev.off()

# Fig3d_堆叠小提琴图（Stacked Violin Plot） 展示特定基因在不同细胞簇中的表达分布差异。--------------------------------------
gene_name <- c(
  "COL6A2",
  "COL4A1",
  "COL14A1",
  "COL5A1",
  "COL5A2",
  "COL8A1"
)

# 执行指定路径下的 R 脚本，加载其中定义的 StackedVlnPlot 函数（用于生成堆叠小提琴图）。
source("2_src/编程猫_src/stackvlion.R")

# 生成图像并处理警告
p <- StackedVlnPlot(
  scRNA, # Seurat对象，包含单细胞数据
  features = gene_name, # 指定要展示的基因（前面定义的6个基因）
  idents = c(5, 1, 4, 3, 0) # 指定要展示的细胞簇及其顺序
) +
  ggsci::scale_fill_lancet() + # 设置填充色为《柳叶刀》期刊配色方案
  scale_y_continuous(name = "Expression Level") # 设置y轴标签为"Expression Level"

# 保存图像
ggsave(
  filename = "Fig/Fig3d_堆叠violion显示指定gene在不同簇中的表达.pdf", # 保存路径和文件名
  plot = p, # 指定要保存的图像对象
  height = 6, width = 5 # 设置图像高度和宽度（单位：英寸）
)
