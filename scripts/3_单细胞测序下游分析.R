# ==============================================
# 单细胞RNA测序下游分析流程
# 版本: 1.0
# 最后更新: 2025-07-12
# ==============================================

rm(list = ls()) # 清空工作空间
# 加载必要的包
pacman::p_load(Seurat, dplyr, data.table, ggplot2, cowplot, patchwork, ggrepel, SingleR, celldex, scater, scran, future, magrittr)
# # 安装SingleR（细胞类型注释核心包）
# BiocManager::install("SingleR", update = FALSE)
#
# # 安装celldex（参考数据集包，SingleR依赖）
# BiocManager::install("celldex", update = FALSE)
#
# # 安装其他可能缺失的Bioconductor包（如scater、scran）
# BiocManager::install(c("scater", "scran"), update = FALSE)

# # 尝试加载celldex包
# library(celldex)
#
# # 检查是否能正常调用包内函数（如获取参考数据集）
# ref <- HumanPrimaryCellAtlasData()
# print(head(ref))  # 若能输出参考数据集信息，说明安装成功

# 设置并行计算
plan("multisession", workers = max(1, parallel::detectCores() - 2))
options(future.globals.maxSize = 8 * 1024^3) # 允许8GB全局变量

# 创建日志函数
log_info <- function(msg) {
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] INFO: ", msg, "\n", sep = "")
}

# 记录分析开始时间
start_time <- Sys.time()
log_info("单细胞RNA测序下游分析开始")

# 0. 加载预处理后的Seurat对象 ----------------------------------------------------

log_info("加载预处理后的Seurat对象...")
scRNA_hq <- readRDS("scRNA_hq.rds") # 请确保已保存此文件
log_info(paste("数据基本情况: 基因数 =", nrow(scRNA_hq), ", 细胞数 =", ncol(scRNA_hq)))

# 1. 原始counts(UMI)计数标准化NormalizeData与特征选择(高变异基因)FindVariableFeatures ------------------------

log_info("开始数据标准化与特征选择...")

# 数据标准化：消除测序深度差异，使不同细胞的基因表达量具有可比性
scRNA_hq <- NormalizeData(
  object = scRNA_hq, # 输入的Seurat对象，经过质量控制的高质量细胞数据
  normalization.method = "LogNormalize", # 标准化方法：对数归一化（最常用）
  scale.factor = 10000 # 缩放因子：将每个细胞的总表达量标准化到该数值
)

# 识别高变异基因HVGs：筛选在细胞间表达差异大的基因，这些基因包含更多生物学信息
scRNA_hq <- FindVariableFeatures(
  object = scRNA_hq, # 输入的Seurat对象（已标准化）
  selection.method = "vst", # 变异基因选择方法：方差稳定变换（推荐）
  nfeatures = 2000 # 保留的高变异基因数量：取前2000个（可调整）
)

# 可视化高变异基因
top10 <- head(VariableFeatures(scRNA_hq), 10)
plot1 <- VariableFeaturePlot(scRNA_hq)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("variable_features.png", plot2, width = 10, height = 8, bg = "white")
log_info("高变异基因分布图已保存至variable_features.png")

# 2. 基因表达量缩放ScaleData与PCA降维RunPCA --------------------------------------------

log_info("开始数据缩放与PCA降维...")

# 数据缩放：消除基因表达量的均值和方差差异，使每个基因对后续分析的影响更均衡
scRNA_hq <- ScaleData(
  object = scRNA_hq, # 输入的Seurat对象（已标准化并筛选高变异基因）
  features = VariableFeatures(scRNA_hq) # 仅对高变异基因进行缩放（节省计算资源）
)

# 主成分分析（PCA）：将高维基因表达数据降维到低维空间，保留主要生物学变异
scRNA_hq <- RunPCA(
  object = scRNA_hq, # 输入的Seurat对象（已缩放）
  features = VariableFeatures(scRNA_hq), # 基于高变异基因进行PCA（这些基因是变异的主要来源）
  npcs = 50 # 计算前50个主成分（PCs）
)

# 3. 确定使用的PC数量（通过Elbow Plot 拐点） ----------------------------------------------
# 生成Elbow Plot并强制修改所有字体属性
p <- ElbowPlot(scRNA_hq, ndims = 50) +
  # 覆盖默认主题，使用ggplot2的基础主题作为起点
  theme_bw() +
  scale_size(range = c(0.5, 1)) +
  theme(
    # 标题（若ElbowPlot显示标题，通常为"Elbow Plot"）
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    # X轴标签（如"PCs"）
    axis.title.x = element_text(size = 16, color = "black"),
    # Y轴标签（如"Standard deviation"）
    axis.title.y = element_text(size = 16, color = "black"),
    # X轴刻度文字（如1,2,...,50）
    axis.text.x = element_text(size = 14, color = "black"),
    # Y轴刻度文字（方差值，如0,5,10...）
    axis.text.y = element_text(size = 14, color = "black"),
    # 确保图例（若存在）的字体大小
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.3)
  )
# 保存图片（确保分辨率足够，避免缩放导致字体变小）
ggsave(
  "PCA_elbow_plot.tiff",
  p,
  width = 6,
  height = 4,
  bg = "white",
  dpi = 300, # 增加分辨率，避免图片压缩导致文字模糊或显小
  compression = "lzw" # 使用LZW压缩，保持图片质量
)
log_info("PCA Elbow图已保存至PCA_elbow_plot.tiff")

# 根据Elbow图结果，选择前20个PC（可根据实际情况调整）
pc_use <- 20
log_info(paste("基于Elbow图，选择前", pc_use, "个PC进行后续分析"))

# 4. 构建细胞邻居图FindNeighbors, 细胞聚类FindClusters与UMAP降维RunUMAP -------------------------------------------------------

log_info("开始细胞聚类与UMAP可视化...")

# 4.1 构建细胞邻居图：基于PCA结果计算细胞间的相似度，为后续聚类做准备
scRNA_hq <- FindNeighbors(
  object = scRNA_hq, # 输入的Seurat对象（已完成PCA）
  dims = 1:pc_use # 使用的主成分范围（如1:15，即前15个PC，pc_use是根据Elbow Plot确定的最佳PC数量）
)

# 4.2 细胞聚类：基于邻居图将相似的细胞归为一类（每个类代表一种潜在的细胞类型）
scRNA_hq <- FindClusters(
  object = scRNA_hq, # 输入已构建邻居图的Seurat对象
  resolution = 0.6 # 聚类分辨率：控制聚类数量（值越大，聚类越细，数量越多）。推荐范围：0.4-1.2
)

# 4.3 UMAP降维：将高维PCA结果进一步压缩到2维空间（UMAP_1和 UMAP_2），使细胞聚类结果可通过散点图直观展示（同类细胞在图中形成紧密的 “簇”）。
# UMAP（Uniform Manifold Approximation and Projection）是一种非线性降维算法，与 PCA 的 “保留全局结构” 不同，它更擅长保留局部邻居关系（相似的细胞在 2 维空间中仍聚集在一起，不同聚类则明显分离）。
scRNA_hq <- RunUMAP(
  object = scRNA_hq, # 输入已完成PCA和聚类的Seurat对象
  dims = 1:pc_use # 使用的主成分范围（与FindNeighbors保持一致，确保基于相同信号）
)

# 4.4 可视化UMAP聚类结果：用散点图展示细胞在2维空间中的分布，不同颜色代表不同聚类
umap_plot <- DimPlot(
  object = scRNA_hq, # 输入已完成UMAP降维的Seurat对象
  reduction = "umap", # 指定可视化的降维结果（这里用UMAP）
  group.by = "seurat_clusters", # 按聚类ID分组（不同聚类用不同颜色）
  label = TRUE, # 在聚类中心显示聚类ID标签（如0、1、2）
  label.size = 5, # 标签字体大小（5为适中值，可根据图片尺寸调整）
  repel = TRUE # 避免标签重叠（自动调整标签位置）
) + ggtitle("UMAP可视化（按聚类）") # 添加图表标题

ggsave("UMAP_clusters.png", umap_plot, width = 6, height = 6, bg = "white")
log_info(paste("聚类完成，共识别出", length(unique(scRNA_hq$seurat_clusters)), "个聚类"))
log_info("UMAP聚类图已保存至UMAP_clusters.png")

# 4.5 tSNE降维（作为UMAP的补充，用于对比不同降维算法的结果）
# tSNE (t-Distributed Stochastic Neighbor Embedding) 是另一种常用的非线性降维方法，
# 特别擅长保留数据的局部结构，但可能扭曲全局距离关系
scRNA_hq <- RunTSNE(
  object = scRNA_hq, # 输入已完成PCA的Seurat对象
  dims = 1:pc_use, # 使用与PCA和UMAP相同的主成分，确保分析一致性
  perplexity = 30, # 控制tSNE对局部和全局结构的平衡，默认值30适用于大多数情况
  # 较小值（如5-10）强调局部细节，较大值（如50-100）更关注全局分布
  # 建议在30-50之间尝试不同值，观察聚类分离效果
  check_duplicates = FALSE # 当数据量较大时，关闭重复检查以提高速度
)

# 可视化tSNE结果
tsne_plot <- DimPlot(
  object = scRNA_hq, # 输入已完成tSNE降维的Seurat对象
  reduction = "tsne", # 指定使用tSNE降维结果进行可视化
  group.by = "seurat_clusters", # 按之前FindClusters()识别的聚类分组着色
  label = TRUE, # 在聚类中心添加聚类标签（如0,1,2...）
  repel = TRUE, # 避免标签重叠（自动调整标签位置）
  pt.size = 0.5 # 控制点的大小（可根据细胞数量调整，细胞越多点越小）
) + ggtitle("tSNE可视化（按聚类）") + # 添加图表标题
  theme(plot.title = element_text(hjust = 0.5)) # 居中对齐标题

# 保存tSNE聚类图为PNG文件
ggsave(
  filename = "tSNE_clusters.png", # 输出文件名
  plot = tsne_plot, # 要保存的ggplot对象
  width = 6, # 图片宽度（英寸）
  height = 6, # 图片高度（英寸）
  dpi = 300, # 分辨率（300dpi适用于大多数出版物）
  bg = "white" # 背景颜色设为白色
)

# 输出分析日志
log_info(paste("tSNE降维完成，结果已保存至tSNE_clusters.png"))
log_info("提示：对比UMAP和tSNE结果可验证聚类稳定性")
log_info(" - 若两者聚类模式相似，说明结果较为可靠")
log_info(" - 若差异较大，可能需要调整tSNE参数或重新评估聚类分辨率")

# 5. 识别每个聚类中高表达的特征性基因FindAllMarkers，差异表达分析与细胞类型注释 --------------------------------------------------------

log_info("开始差异表达分析与细胞类型注释...")

# 5.1 寻找每个聚类的标记基因：通过差异表达分析，识别每个聚类中高表达的特征性基因
cluster_markers <- FindAllMarkers(
  object = scRNA_hq, # 输入已完成聚类的Seurat对象
  only.pos = TRUE, # 仅保留在聚类中相对其他聚类高表达的阳性标记（排除低表达基因）
  min.pct = 0.25, # 要求基因至少在25%的目标聚类细胞中表达（过滤低丰度基因）
  logfc.threshold = 0.25 # 要求基因在目标聚类中的表达量比其他聚类至少高log2(1.25)≈0.22
)

top10 <- cluster_markers %>%
  group_by(cluster) %>% # 按聚类分组
  top_n(n = 10, wt = avg_log2FC) # 每个聚类取前10个基因（按avg_log2FC排序）

# 可视化每个聚类的top10标记基因
marker_plot <- DotPlot(
  object = scRNA_hq, # 输入已完成聚类的Seurat对象
  features = top10$gene, # 使用前10个标记基因
  group.by = "seurat_clusters", # 按聚类分组显示
  dot.scale = 8, # 点的大小（表示表达量）
  col.min = 0, # 最小颜色值（表示表达量）
  col.max = 3 # 最大颜色值（表示表达量）
) + ggtitle("Top 10 Markers per Cluster") + # 添加标题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # X轴标签倾斜45度，便于阅读

marker_plot
ggsave("top10_markers_per_cluster.png", marker_plot, width = 10, height = 6, bg = "white")


# 绘制DoHeatmap图 ------------------------------------------------------------
# 1. 获取所有有效基因名称（Seurat v5兼容方式）
Assays(scRNA_hq)
# 获取所有assay名称（返回字符向量）
assay_names <- names(scRNA_hq)
print(assay_names) # ["RNA"     "RNA_nn"  "RNA_snn" "pca"     "umap"    "tsne"] ，"RNA"：存储基因表达数据。"RNA_nn" 和 "RNA_snn"：存储近邻图（nearest neighbor graph）和共享近邻图（shared nearest neighbor graph）。"pca", "umap", "tsne"：存储降维结果（如 PCA、UMAP、t-SNE）。

all_genes <- rownames(scRNA_hq[["RNA"]]) # [["RNA"]] 表示提取名为"RNA"的assay（分析组），其中存储了基因表达的原始计数或标准化数据。旧版是rownames(scRNA_hq@assays$RNA@counts)

# 2. 验证top10基因是否都存在于数据中
valid_features <- top10$gene[top10$gene %in% all_genes]

# 3. 检查是否有基因被过滤掉
if (length(valid_features) < length(top10$gene)) {
  missing_genes <- setdiff(top10$gene, valid_features)
  message("以下基因不存在于表达矩阵中，已被过滤: ", paste(missing_genes, collapse = ", "))
}

# 4. 确保数据已正确缩放（使用valid_features）
scRNA_hq <- ScaleData(scRNA_hq, features = valid_features)

# 5. 使用验证后的有效基因列表绘制热图
p <- DoHeatmap(
  object = scRNA_hq,
  features = valid_features, # 使用过滤后的有效基因
  group.by = "seurat_clusters",
  size = 3,
  draw.lines = TRUE
) + ggtitle("Top 10 Valid Markers Heatmap")

# 6. 保存热图
ggsave("top10_valid_markers_heatmap.png", plot = p, width = 10, height = 6, bg = "white")


# 修改热图的颜色scale_fill_gradientn() -----------------------------------------------------------------
library(ggplot2)

p1 <- DoHeatmap(
  object = scRNA_hq,
  features = valid_features, # 使用过滤后的有效基因
  group.by = "seurat_clusters",
  size = 3,
  draw.lines = TRUE
) + scale_fill_gradientn(
  colors = c("blue", "white", "red"), # 从蓝色到白色再到红色的渐变
  values = scales::rescale(c(-2, 0, 2)), # 确保颜色映射到-2到2的范围
  limits = c(-2, 2) # 设置颜色映射的范围
) + ggtitle("Top 10 Valid Markers Heatmap with Custom Colors")

p1

# 可视化特定基因在不同细胞聚类中的表达情况（如免疫细胞标记物）FeaturePlot与VlnPlot --------------------------------------------

# 定义一组感兴趣的基因（如免疫细胞标记物）CD4+ T 细胞：高表达 CD4、FOXP3，CD8+ T 细胞：高表达 CD8A、CD8B、CD3D/E，B 细胞：高表达 CD19、MS4A1（CD20），单核细胞 / 巨噬细胞：高表达 CD14、LYZ。
select_genes <- c("CD4", "FOXP3", "CD8A", "CD8B", "CD3D", "CD3E", "CD19", "MS4A1", "CD14", "LYZ")

# FeaturePlot：绘制特定基因在UMAP空间中的表达分布，颜色越深表示基因表达水平越高
p2 <- FeaturePlot(
  object = scRNA_hq,
  features = select_genes, # 使用选择的基因
  reduction = "umap", # reduction指定使用哪种降维结果（如 UMAP、t-SNE、PCA），使用UMAP降维结果展示每个基因在细胞中的表达分布
  label = TRUE, # 在图中标记细胞类型
  ncol = 3, # 每行显示3个基因
) + ggtitle("Selected Genes Heatmap with Custom Colors")

p2

# VlnPlot：绘制小提琴图展示基因在不同聚类中的表达
p3 <- VlnPlot(
  object = scRNA_hq, # Seurat对象
  features = select_genes, # 使用选择的基因
  group.by = "seurat_clusters", # 按聚类分组（展示每个基因在不同聚类中的表达）
  pt.size = 0.1, # 点的大小（控制每个数据点的显示大小）
  combine = TRUE # # 将所有基因的小提琴图合并为一个图
) + ggtitle("Selected Genes Violin Plot") # 添加图标题

p3


# 保存标记基因
write.csv(cluster_markers, "cluster_markers.csv", row.names = FALSE)
# CSV文件关键列解析：
# gene: 基因名称（如CD4、FOXP3）
# cluster: 该基因所属的聚类（如0、1、2...）
# p_val: 差异表达的p值（越小越显著）
# avg_log2FC: 对数转换后的表达量差异倍数（正值表示在该聚类中高表达）
# pct.1: 基因在目标聚类中的表达比例（如0.8表示80%的细胞表达该基因）
# pct.2: 基因在其他聚类中的表达比例
# p_val_adj: 校正后的p值（多重检验校正，通常要求<0.05）
log_info("聚类标记基因已保存至cluster_markers.csv")

# 使用dplyr包按聚类分组，提取每个聚类的top5标记基因（按表达量差异降序排列）
top_markers <- cluster_markers %>%
  group_by(cluster) %>% # 按聚类分组
  arrange(desc(avg_log2FC)) %>% # 每个聚类内按avg_log2FC降序排列
  slice_head(n = 5) # 取前5个基因

print(top_markers)


# 6. 细胞类型注释（SingleR） ---------------------------------------------------------

log_info("开始细胞类型注释...")

###### 方法1：基于已知标记基因的手动注释（示例映射，需根据实际数据调整）#####
# cell_type_mapping <- c(
#   "0" = "CD4+ T细胞",
#   "1" = "CD8+ T细胞",
#   "2" = "B细胞",
#   "3" = "单核细胞",
#   "4" = "自然杀伤细胞",
#   "5" = "树突状细胞",
#   "6" = "巨噬细胞",
#   "7" = "内皮细胞",
#   "8" = "成纤维细胞",
#   "9" = "浆细胞"
#   # 根据实际标记基因添加更多映射...
# )
#
# # 应用细胞类型注释
# scRNA_hq$manual_cell_type <- cell_type_mapping[as.character(scRNA_hq$seurat_clusters)]
# scRNA_hq$manual_cell_type <- factor(scRNA_hq$manual_cell_type)

##### 方法2：使用SingleR的自动注释 #####
log_info("执行SingleR自动细胞类型注释...")
ref <- HumanPrimaryCellAtlasData() # 加载人类原代细胞参考数据集，是一个SummarizedExperiment对象（ Bioconductor 中用于存储高维生物数据（如 RNA-seq、单细胞数据）的标准容器）。SummarizedExperiment 对象包含5个核心组件：colData：样本元数据（如细胞类型、聚类信息）。rowData 或 NAMES：行（基因名）。assays：表达矩阵。elementMetadata: 特征元数据（通常为空）。metadata: 数据集整体元数据。

str(ref)

# 获取样本元数据（细胞注释信息）
cell_metadata <- colData(ref)

# 查看基本信息
print(cell_metadata) # 显示基本结构
dim(cell_metadata) # 查看行数（细胞数）和列数（元数据字段数）
names(cell_metadata) # 查看元数据包含哪些字段，"label.main"宏观分类-主要细胞类型（如 "DC"、"Monocyte"），是简化的分类。 "label.fine"精细细胞类型（如 "DC:monocyte-derived:immature"），更详细的亚群描述。"label.ont" 标准化的细胞类型标识符-细胞类型的 ontology ID（如 "CL:0000840"），关联到标准化的细胞类型术语（如 Cell Ontology细胞本体论）。

# 查看具体内容
head(cell_metadata) # 显示前几个细胞的元数据
table(cell_metadata$label.main) # 统计每种主要细胞类型的数量

# 准备SingleR输入：提取Seurat对象中的原始计数矩阵（slot = "counts"）
test <- GetAssayData(scRNA_hq, slot = "counts")
class(test) # 确认test是一个矩阵（matrix）或稀疏矩阵（dgCMatrix）
dim(test) # 查看测试数据的维度（基因数 x 细胞数）
rownames(test)[1:10] # 查看基因名称（行名）
colnames(test)[1:10] # 查看细胞名称（列名）
subset_counts <- test[1:5, 1:5] # 提取前100个基因和前10个细胞进行测试
subset_counts_dense <- as.data.frame(as.matrix(subset_counts)) # 将稀疏矩阵转换为数据框
print(subset_counts) # 打印前5个基因和前5个细胞的计数矩阵

## 执行SingleR注释，将测试数据中的每个细胞分类到参考数据的已知细胞类型
singleR_pred <- SingleR(
  test = test, # 待注释的测试数据（原始计数矩阵）
  ref = ref, # 参考数据集
  labels = ref$label.main # 参考数据中的细胞类型标签（主要分类）
)

## 将注释结果添加到Seurat对象
scRNA_hq$singleR_cell_type <- singleR_pred$labels

# 查看元数据的前几行
head(scRNA_hq@meta.data)

# 获取所有唯一的细胞类型
cell_types <- unique(scRNA_hq$singleR_cell_type)
print(cell_types) # 打印所有唯一的细胞类型

# 设置Seurat对象的细胞类型标识符（注意：Idents首字母大写）
Idents(scRNA_hq) <- scRNA_hq$singleR_cell_type
# 验证设置是否成功
head(Idents(scRNA_hq)) # 输出当前的细胞身份
table(Idents(scRNA_hq)) # 统计每种细胞类型的数量

# 可视化两种注释结果
# manual_plot <- DimPlot(scRNA_hq, reduction = "umap", group.by = "manual_cell_type", label = TRUE, repel = TRUE) + ggtitle("手动注释")
auto_plot <- DimPlot(scRNA_hq, reduction = "umap", group.by = "singleR_cell_type", label = TRUE, repel = TRUE) + ggtitle("SingleR自动注释")

auto_plot

ggsave("SingleR_annotation.tiff", auto_plot, width = 6, height = 6, bg = "white", dpi = 300, compression = "lzw")

# annotation_comparison <- manual_plot | auto_plot
# ggsave("cell_type_annotation_comparison.png", annotation_comparison, width = 18, height = 8)
log_info("细胞类型注释完成")
log_info("细胞类型注释对比图已保存至cell_type_annotation_comparison.png")

# ---------------------
# 6. 保存最终结果
# ---------------------
saveRDS(scRNA_hq, file = "scRNA_hq_annotated.rds")
log_info("完整分析结果已保存至scRNA_hq_annotated.rds")

# 记录分析结束时间
end_time <- Sys.time()
log_info(paste("分析完成，总耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))
