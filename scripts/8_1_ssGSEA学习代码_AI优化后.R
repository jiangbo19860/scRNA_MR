#-----------------------------------------------------------------------------
# 单细胞RNA-seq数据的通路活性分析与差异检验流程
# 功能：使用ssGSEA算法计算Hallmark通路活性，并分析肿瘤与非恶性细胞的差异
#-----------------------------------------------------------------------------

# 清空环境变量
rm(list = ls())

# 加载必要的包（自动安装缺失包）
pacman::p_load(
  here, # 路径管理
  Seurat, # 单细胞数据分析
  ggplot2, # 绘图
  clusterProfiler, # 富集分析
  GSVA, # 基因集变异分析
  tidyverse, # 数据处理
  BiocParallel # 并行计算
)

# 设置工作目录
here()

#-----------------------------------------------------------------------------
# 1. 数据读取与预处理
#-----------------------------------------------------------------------------
# 读取Hallmark基因集（GMT格式）
hallmark <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt"))
# 去除通路名称前缀"HALLMARK_"
hallmark$term <- gsub("HALLMARK_", "", hallmark$term)

# 将数据框转换为通路-基因列表格式
hallmark.list <- hallmark %>%
  split(.$term) %>% # 按通路名称分组，split() 函数的作用是把数据框或者向量按照某个因子进行分组，最终返回一个列表。split(x, f). x：输入对象，一般是数据框或者向量。f：一个因子或者可转换为因子的对象，它定义了分组的依据。
  lapply("[[", 2) # 从每个分组后的子数据框中提取第二列（基因名称列），最终将通路-基因的对应关系从 "数据框列表"（每个通路对应一个子数据框，包含2列信息（通路名和基因名） 转换为 "纯基因列表"（每个通路直接对应一个基因名称向量（仅保留核心信息））。

# 加载单细胞RNA-seq数据（Seurat对象）
load(here("1_data", "E-MTAB-6149", "scRNA_endo.Rds"))
# 提取标准化后的表达矩阵（log2(counts+1)）
expr_norm <- as.matrix(scRNA@assays$RNA@data)

#-----------------------------------------------------------------------------
# 2. ssGSEA通路活性分析
#-----------------------------------------------------------------------------
# 过滤常量基因（关键步骤：避免因全零或恒定表达值基因导致的错误）
is_constant <- apply(expr_norm, 1, function(x) length(unique(x)) == 1)
expr_filtered <- expr_norm[!is_constant, ]

# 查看过滤效果
cat("过滤前基因数:", nrow(expr_norm), "\n")
cat("过滤后基因数:", nrow(expr_filtered), "\n")
cat("被过滤的常量基因数:", sum(is_constant), "\n")

# 检查表达矩阵与基因集的基因名称匹配情况（关键步骤：避免基因映射失败）
matrix_genes <- rownames(expr_filtered)
all_geneset_genes <- unique(unlist(hallmark.list))

# 计算交集比例（越高越好，理想情况>80%）
common_genes <- intersect(matrix_genes, all_geneset_genes)
cat("表达矩阵与基因集的共同基因数:", length(common_genes), "\n")
cat("基因集可映射比例:", length(common_genes) / length(all_geneset_genes) * 100, "%\n")

# 清理基因集：移除在表达矩阵中没有匹配基因的通路（关键步骤：避免空基因集）
hallmark_clean <- lapply(hallmark.list, function(genes) {
  intersect(genes, matrix_genes) # 只保留在表达矩阵中存在的基因
})

# 删除空基因集（如果某个通路的所有基因都不在表达矩阵中）
hallmark_clean <- hallmark_clean[sapply(hallmark_clean, length) > 0]
cat("清理后保留的通路数:", length(hallmark_clean), "（原始通路数:", length(hallmark.list), "）\n")

# 构建ssGSEA参数对象（修正：使用过滤后的表达矩阵和清理后的基因集）
ssgsea_param <- ssgseaParam(
  expr = expr_filtered, # 使用过滤后的表达矩阵
  geneSets = hallmark_clean # 使用清理后的基因集
)

# 设置并行计算（根据CPU核心数调整）
register(MulticoreParam(workers = 3))

# 执行ssGSEA分析（通路活性计算）
# 注意：结果为通路×细胞的矩阵（行=通路，列=细胞）
gsva_result <- gsva(
  param = ssgsea_param, # 分析参数
  BPPARAM = bpparam(), # 并行计算参数
  verbose = FALSE # 静默运行
)

# 检查结果维度
cat("GSVA结果维度（通路数×细胞数）:", dim(gsva_result), "\n")

# 将结果整合到Seurat对象的meta.data中（确保转置为细胞×通路）
scRNA@meta.data <- cbind(scRNA@meta.data, t(gsva_result))

# 去重列（防止重复列名导致后续分析错误）
duplicated_cols <- duplicated(colnames(scRNA@meta.data))
scRNA@meta.data <- scRNA@meta.data[, !duplicated_cols]

# 验证结果
cat("成功将", nrow(gsva_result), "个通路的活性分数添加到Seurat对象中\n")


#  3. 差异通路分析（Wilcoxon秩和检验）-------------------------------
# 获取所有通路名称
pathways <- rownames(gsva_result)

# 提取用于差异分析的数据（分组信息+通路活性）
meta_data <- scRNA@meta.data %>%
  select(all_of(c("CellFromTumor", pathways)))

# 定义Wilcoxon检验函数（对单个通路进行分析）
test_wilcox <- function(pathway) {
  # 提取特定通路的活性数据和分组信息
  df <- meta_data %>%
    select(CellFromTumor, all_of(pathway)) %>% # 选择列
    drop_na() %>% # 去除缺失值
    mutate(
      # 将分组变量转换为因子并设置标签
      CellFromTumor = factor(
        CellFromTumor,
        levels = c("0", "1"), # 原始值
        labels = c("Non-malignant", "Tumor") # 新标签
      )
    )

  # 执行Wilcoxon秩和检验
  wilcox_result <- wilcox.test(
    as.formula(paste(pathway, "~ CellFromTumor")), # 公式：通路活性 ~ 分组
    data = df
  )

  # 返回结果数据框
  return(data.frame(
    pathway = pathway,
    statistic = wilcox_result$statistic,
    p_value = wilcox_result$p.value,
    stringsAsFactors = FALSE
  ))
}

# 并行执行所有通路的差异检验
register(MulticoreParam(workers = 4))
wilcox_results <- bplapply(
  X = pathways, # 通路列表
  FUN = test_wilcox, # 检验函数
  BPPARAM = bpparam() # 并行参数
) %>%
  bind_rows() %>% # 合并结果
  mutate(
    fdr = p.adjust(p_value, method = "fdr") # FDR多重检验校正
  ) %>%
  arrange(fdr) # 按FDR排序

# 筛选显著差异通路（FDR < 0.05）
wilcox_sig <- wilcox_results %>%
  filter(fdr < 0.05) %>% # 显著性阈值
  arrange(desc(statistic)) # 按统计量降序排列

#-----------------------------------------------------------------------------
# 4. 结果可视化与保存
#-----------------------------------------------------------------------------
# 设置输出目录
output_dir <- here("3_outputs", "ssGSEA") # 用逗号分隔路径组件，无需手动加/

# 创建目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE) # 自动创建多级目录
  cat("成功创建子文件夹：", output_dir, "\n")
} else {
  cat("子文件夹已存在：", output_dir, "\n")
}

#-----------------------------------------------------------------------------
# 4.1 保存表格结果
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 定义数值格式化函数：小数点后3位至少有一个>0，则保留3位小数，否则使用科学计数法保留2位有效数字
#-----------------------------------------------------------------------------
format_num <- function(x) {
  if (!is.numeric(x)) {
    return(as.character(x))
  }

  # 将数值转换为字符，保留足够的小数位
  x_char <- format(x, nsmall = 10)

  # 提取小数点后前3位（如果有小数点）
  decimal_part <- sub(".*\\.", "", x_char)
  decimal_part <- substr(decimal_part, 1, 3)

  # 判断条件：小数点后3位中至少有一个大于0
  has_nonzero <- any(as.numeric(strsplit(decimal_part, "")[[1]]) > 0)

  # 根据条件选择格式
  if (has_nonzero) {
    return(format(round(x, 3), nsmall = 3)) # 保留3位小数
  } else {
    return(format(x, scientific = TRUE, digits = 2)) # 科学计数法，2位有效数字
  }
}

#-----------------------------------------------------------------------------
# 保存分析结果（应用格式化函数）
#-----------------------------------------------------------------------------
# 保存显著差异通路结果
wilcox_sig_formatted <- data.frame(lapply(wilcox_sig, format_num))
write.csv(
  wilcox_sig_formatted,
  file.path(output_dir, "significant_pathways.csv"),
  row.names = FALSE
)

# 保存所有通路的检验结果
wilcox_results_formatted <- data.frame(lapply(wilcox_results, format_num))
write.csv(
  wilcox_results_formatted,
  file.path(output_dir, "all_pathways_wilcox_results.csv"),
  row.names = FALSE
)

# 保存通路活性矩阵（注意：row.names=TRUE保留通路名称）
gsva_result_formatted <- data.frame(lapply(gsva_result, format_num))
write.csv(
  gsva_result_formatted,
  file.path(output_dir, "gsva_pathway_activities.csv"),
  row.names = TRUE
)

#-----------------------------------------------------------------------------
# 可选：保存其他结果（如果有）
#-----------------------------------------------------------------------------
# 保存通路富集热图数据
if (exists("pathway_heatmap_data")) {
  pathway_heatmap_data_formatted <- data.frame(lapply(pathway_heatmap_data, format_num))
  write.csv(
    pathway_heatmap_data_formatted,
    file.path(output_dir, "pathway_heatmap_data.csv"),
    row.names = TRUE
  )
}

#-----------------------------------------------------------------------------
# 4.2 保存可视化结果
#-----------------------------------------------------------------------------
# 4.2.1 小提琴图：展示显著通路在两组间的活性分布
if (nrow(wilcox_sig) > 0) {
  # 选择前5个最显著的通路
  top_pathways <- head(wilcox_sig$pathway, 5)

  # 构建小提琴图
  p_vln <- VlnPlot(
    scRNA,
    features = top_pathways,
    group.by = "CellFromTumor", # 按分组展示
    pt.size = 0.2, # 点大小
    cols = c("royalblue", "firebrick"), # 颜色方案
    ncol = 1 # 按列排列
  ) +
    theme(
      axis.text.x = element_text(size = 10), # x轴文本大小
      plot.title = element_text(hjust = 0.5) # 标题居中
    )

  # 保存小提琴图
  ggsave(
    file.path(output_dir, "pathway_activity_vlnplot.pdf"),
    plot = p_vln,
    width = 6,
    height = 4 * length(top_pathways) # 根据通路数量调整高度
  )
  print(p_vln) # 显示小提琴图
}

# 4.2.2 火山图：展示所有通路的差异显著性
p_volcano <- ggplot(wilcox_results, aes(x = statistic, y = -log10(fdr))) +
  geom_point(
    aes(color = fdr < 0.05), # 颜色根据显著性区分
    size = 2, # 点大小
    alpha = 0.7 # 透明度
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") + # 垂直参考线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") + # 水平参考线
  scale_color_manual(values = c("gray50", "red")) + # 颜色方案
  labs(
    x = "Wilcoxon Statistic (Tumor vs Non-malignant)",
    y = "-log10(FDR)",
    title = "Pathway Activity Differences"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # 不显示图例
print(p_volcano)

# 保存火山图
ggsave(
  file.path(output_dir, "pathway_activity_volcano.pdf"),
  plot = p_volcano,
  width = 7,
  height = 5
)

# 4.2.3 UMAP特征图：展示特定通路在UMAP上的分布
p_umap <- FeaturePlot(scRNA, features = "ADIPOGENESIS", reduction = "umap")

# 保存UMAP图
ggsave(
  file.path(output_dir, "ADIPOGENESIS_umap_featureplot.pdf"),
  plot = p_umap,
  width = 6,
  height = 5
)

#-----------------------------------------------------------------------------
# 5. 分析完成提示
#-----------------------------------------------------------------------------
cat("分析完成！结果已保存至:", output_dir, "\n")
