### 单细胞RNA-seq数据的基因集富集分析ssGSEA，用ssgsea算法分析 hallmark 基因集在单细胞中的通路活性。
# 输入：标准化的单细胞表达矩阵（expr_norm）和 Hallmark 基因集（hallmark.list）；
# 输出：每个Hallmark通路在每个细胞中的富集分数（活性值）。

rm(list = ls())

pacman::p_load(here, Seurat, ggplot2, clusterProfiler, GSVA, tidyverse, BiocParallel)
here()

hallmark <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt"))
colnames(hallmark)
head(hallmark)
str(hallmark)

hallmark$term <- gsub("HALLMARK_", "", hallmark$term)

hallmark.list <- hallmark %>%
  split(.$term) %>% # 把数据框按照term列进行分组，生成一个子数据框列表list。.代表的是前一个表达式所产生的结果.$是从数据框或者列表中提取特定列。
  lapply("[[", 2) # 把"[["函数作为参数传递给lapply，同时将2作为额外参数。这就相当于对列表中的每个子数据框都执行"[["(x, 2)操作，也就是提取每个子数据框的第二列（gene symbol），生成一个列表。lapply是对列表中的每个元素执行相同的操作，并且最终返回一个列表。


load(here("1_data", "E-MTAB-6149", "scRNA_endo.Rds"))
scRNA
DimPlot(scRNA)
DefaultDimReduc(scRNA) # 输出类似 "umap" 或 "tsne"
meta <- scRNA@meta.data
head(meta)
expr <- as.matrix(scRNA@assays$RNA@counts)
dim(expr)

# 使用Seurat标准化后的表达矩阵（log归一化，避免counts的量纲影响）
expr_norm <- as.matrix(scRNA@assays$RNA@data) # Seurat的data层为log2(counts+1)归一化后的数据

# 查看Seurat对象的命令历史（记录了执行过的主要函数）
# scRNA@commands  # 查看Seurat对象的命令历史（记录了执行过的主要函数）

# 检查是否已经进行过log归一化 ------------------------------------------------------
# 提取data层数据的部分值查看
head(expr_norm[1:5, 1:5]) # 查看前5个基因、前5个细胞的表达值
# 统计数据范围（log2归一化后的值通常在0~10左右，极少超过20）
range(expr_norm)
# 检查是否有负值（log2归一化后理论上不会有负值，因counts≥0，log2(counts+1)≥0）
any(expr_norm < 0)


# 运行ssgsea，将基因水平的单细胞表达数据，转换为通路水平的活性分数，得到每个细胞中 Hallmark 通路的活性分数矩阵。 ----------------------------------------------------------------
# 第一步：通过ssgseaParam()构造函数，将表达数据矩阵和基因集，整合为算法可识别的参数对象，明确使用ssGSEA算法。隐含参数：最新版本的 GSVA 会自动使用默认设置（如kcdf="Gaussian"，适合 log 归一化数据；absRanking=TRUE，基于绝对表达值排序）。
ssgsea_param <- ssgseaParam(
  expr = expr_norm, # 输入：标准化后的表达矩阵（基因×细胞）
  geneSets = hallmark.list # 输入：Hallmark基因集列表（50个核心通路）
)

# 通过BiocParallel包中的MulticoreParam()启用并行计算，加速通路富集分析（尤其适合 1592 个细胞的大规模数据）。
register(MulticoreParam(workers = 3)) # 3个核心，可根据电脑配置调整

# 重要：运行ssGSEA分析, 使用参数对象调用gsva(), 得到的gsva_result是一个通路 × 细胞的矩阵（基因集富集分数矩阵)。行：50 个 Hallmark 通路（如ADIPOGENESIS、APOPTOSIS等）；列：1592 个细胞（与输入表达矩阵的细胞一一对应）；值：通路富集分数（数值越大，该通路在对应细胞中的活性越高）。
gsva_result <- gsva(
  param = ssgsea_param, # 传入参数对象（包含数据和算法设置）
  BPPARAM = bpparam(), # 应用并行计算参数（使用3个核心）
  verbose = FALSE # 不输出中间过程（静默运行）
)
class(gsva_result) # 检查结果类型，应该是matrix，矩阵中所有元素必须是同一数据类型（例如全部是数值型）。数据框的每一列可以是不同的数据类型（例如一列是数值、一列是字符串、一列是因子）。

# 将结果整合到Seurat对象的meta.data中，用于可视化（需转置为细胞×通路）
scRNA@meta.data <- cbind(scRNA@meta.data, t(gsva_result))

# 查看新增的通路列（如"ADIPOGENESIS"）
colnames(scRNA@meta.data)

meta <- scRNA@meta.data
unique(meta$ClusterName)
unique(meta$CellFromTumor)

# 在UMAP上可视化某通路活性（如"ADIPOGENESIS"）
FeaturePlot(scRNA, features = "ADIPOGENESIS", reduction = "umap")


# 基于 Wilcoxon 秩和检验，识别在肿瘤细胞（Tumor）和非恶性细胞（Non-malignant）中显著差异激活的通路，并 --------
# ----------------------------
# 步骤1：去除meta.data中重复的通路列
# ----------------------------
# 查看重复的通路列（找出重复的列名）
duplicated_cols <- duplicated(colnames(scRNA@meta.data))

print(duplicated_cols)

# 保留第一次出现的列，去除重复列
scRNA@meta.data <- scRNA@meta.data[, !duplicated_cols]

# ----------------------------
# 步骤2：提取包含分组和通路活性的数据
# ----------------------------
# 通路名称是gsva_result的行名（已整合到meta.data中）
pathways <- rownames(gsva_result)

# 提取数据：包含CellFromTumor和所有通路列
meta_data <- scRNA@meta.data %>%
  select(all_of(c("CellFromTumor", pathways))) # 确保通路列正确提取

# ----------------------------
# 步骤3：定义Wilcoxon检验函数（含分组标签转换）
# ----------------------------
test_wilcox <- function(pathway) {
  # 提取该通路的活性数据和分组信息
  df <- meta_data %>%
    select(CellFromTumor, all_of(pathway)) %>% # 提取分组和当前通路
    drop_na() %>% # 去除缺失值
    mutate(
      CellFromTumor = factor(
        CellFromTumor,
        levels = c("0", "1"), # 原始值：0=非恶性，1=肿瘤
        labels = c("Non-malignant", "Tumor") # 新标签
      )
    )

  # 执行Wilcoxon秩和检验（两组比较）
  wilcox_result <- wilcox.test(
    as.formula(paste(pathway, "~ CellFromTumor")), # 通路活性 ~ 分组
    data = df
  )

  # 返回结果：通路名称、检验统计量、原始p值
  return(data.frame(
    pathway = pathway,
    statistic = wilcox_result$statistic,
    p_value = wilcox_result$p.value,
    stringsAsFactors = FALSE
  ))
}

# ----------------------------
# 步骤4：批量检验所有通路的组间差异
# ----------------------------
# 并行计算（可选，加速分析）
library(BiocParallel)
register(MulticoreParam(workers = 4)) # 使用4个核心

# 对每个通路执行检验
wilcox_results <- bplapply(
  X = pathways,
  FUN = test_wilcox,
  BPPARAM = bpparam() # 应用并行参数
) %>%
  bind_rows() %>% # 合并结果为数据框
  mutate(
    fdr = p.adjust(p_value, method = "fdr") # FDR校正
  ) %>%
  arrange(fdr) # 按FDR从小到大排序

# ----------------------------
# 步骤5：筛选显著差异通路（FDR < 0.05）
# ----------------------------
wilcox_sig <- wilcox_results %>%
  filter(fdr < 0.05) %>% # 显著阈值
  arrange(desc(statistic)) # 按统计量排序（正→Tumor中活性更高，负→Non-malignant中更高）

# 查看显著通路
print(wilcox_sig)

# ----------------------------
# 步骤6：可视化显著通路（示例）
# ----------------------------
if (nrow(wilcox_sig) > 0) {
  # 选择前5个最显著的通路可视化
  top_pathways <- head(wilcox_sig$pathway, 5)

  # 小提琴图：展示通路活性在两组间的分布
  VlnPlot(
    scRNA,
    features = top_pathways,
    group.by = "CellFromTumor", # 按分组展示
    pt.size = 0.5, # 点大小
    cols = c("royalblue", "firebrick"), # 颜色：Non-malignant=蓝，Tumor=红
    ncol = 1 # 按列排列
  ) +
    theme(
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5) # 标题居中
    )

  # 火山图：展示所有通路的差异显著性
  ggplot(wilcox_results, aes(x = statistic, y = -log10(fdr))) +
    geom_point(
      aes(color = fdr < 0.05),
      size = 2,
      alpha = 0.7
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("gray50", "red")) +
    labs(
      x = "Wilcoxon Statistic (Tumor vs Non-malignant)",
      y = "-log10(FDR)",
      title = "Pathway Activity Differences"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}


# 提取包含通路活性和CellFromTumor的数据
meta_data <- scRNA@meta.data %>%
  select(all_of(c("CellFromTumor", names(gsva_result)))) # 保留分组和通路列
pathways <- rownames(gsva_result)
# 定义修改后的检验函数（含标签替换）
test_wilcox <- function(pathway) {
  df <- meta_data %>%
    select(CellFromTumor, all_of(pathway)) %>%
    drop_na() %>%
    mutate(CellFromTumor = factor(
      CellFromTumor,
      levels = c("0", "1"), # 原始值
      labels = c("Non-malignant", "Tumor") # 新标签
    ))

  wilcox_result <- wilcox.test(as.formula(paste(pathway, "~ CellFromTumor")), data = df)
  return(data.frame(
    pathway = pathway,
    statistic = wilcox_result$statistic,
    p_value = wilcox_result$p.value
  ))
}

# 批量检验（以所有通路为例）
wilcox_results <- lapply(pathways, test_wilcox) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(fdr)

# 查看显著通路
wilcox_sig <- wilcox_results %>% filter(fdr < 0.05)
print(wilcox_sig)
