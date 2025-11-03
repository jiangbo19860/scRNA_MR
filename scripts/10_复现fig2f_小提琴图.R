rm(list = ls())
pacman::p_load(Seurat, dplyr, limma, stringr, GSVA, clusterProfiler, data.table, ggplot2, cowplot, patchwork, ggrepel, hdf5r, BiocParallel)
here()
load(here("3_outputs", "E-MTAB-6149_outputs", "scRNA_endo.Rds"))

# 完整的处理流程
hallmark <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt")) %>%
  mutate(term = gsub("HALLMARK_", "", term)) %>% # 清理通路名称
  split(.$term) %>% # 按通路分组
  lapply("[[", 2) # 转换为通路-基因列表

# 提取标准化后的表达矩阵（Seurat 5.0+版本）
expr_norm <- GetAssayData(scRNA, assay = "RNA", layer = "data") %>%
  as.matrix()
dim(expr_norm) # 查看表达矩阵维度（基因数 x 细胞数）)


# 过滤常量基因
is_constant <- apply(expr_norm, 1, function(x) length(unique(x)) == 1)
expr_filtered <- expr_norm[!is_constant, ]
dim(expr_filtered) # 查看过滤后的表达矩阵前几行
expr_filtered[1:5, 1:5] # 打印前5行5列数据

# 检查表达矩阵与基因集的基因名称匹配情况（关键步骤：避免基因映射失败）
matrix_genes <- rownames(expr_filtered)
all_geneset_genes <- unique(unlist(hallmark))

# 计算交集比例（越高越好，理想情况>80%）
common_genes <- intersect(matrix_genes, all_geneset_genes)
cat("表达矩阵与基因集的共同基因数:", length(common_genes), "\n")
cat("基因集可映射比例:", length(common_genes) / length(all_geneset_genes) * 100, "%\n")

# 清理基因集：移除在表达矩阵中没有匹配基因的通路（关键步骤：避免空基因集）
hallmark_clean <- lapply(hallmark, function(genes) {
  intersect(genes, matrix_genes) # 只保留在表达矩阵中存在的基因
})

# 删除空基因集（如果某个通路的所有基因都不在表达矩阵中）
hallmark_clean <- hallmark_clean[sapply(hallmark_clean, length) > 0]
cat("清理后保留的通路数:", length(hallmark_clean), "（原始通路数:", length(hallmark), "）\n")

# 构建ssGSEA参数对象（使用过滤后的表达矩阵和清理后的基因集）
ssgsea_param <- ssgseaParam(
  expr = expr_filtered, # 使用过滤后的表达矩阵
  geneSets = hallmark_clean # 使用清理后的基因集
)

# 设置并行计算（根据CPU核心数调整）
register(MulticoreParam(workers = 3))

# 执行ssGSEA分析（通路活性计算）
# 注意：结果gsva_result为通路×细胞的矩阵（行=通路，列=细胞）
gsva_result <- gsva(
  param = ssgsea_param, # 分析参数
  BPPARAM = bpparam(), # 并行计算参数
  verbose = FALSE # 静默运行
)

# 检查结果维度
cat("GSVA结果维度（通路数×细胞数）:", dim(gsva_result), "\n")

# 将结果整合到Seurat对象的meta.data中（确保转置为细胞×通路）
scRNA@meta.data <- cbind(scRNA@meta.data, t(gsva_result))
dim(scRNA)
meta <- scRNA@meta.data

# 去重列（防止重复列名导致后续分析错误）
duplicated_cols <- duplicated(colnames(scRNA@meta.data))
scRNA@meta.data <- scRNA@meta.data[, !duplicated_cols]
dim(scRNA)

# 验证结果
cat("成功将", nrow(gsva_result), "个通路的活性分数添加到Seurat对象中\n")


# 绘制多个小提琴图 ----------------------------------------------------------------
gene_name <- c(
  "HSPG2",
  "ANGPT2",
  "HIF1A",
  "MMP2",
  "CTGF",
  "NOTCH1"
)

# 初始化数据框，包含细胞聚类信息
data_tmp <- data.frame(
  group = scRNA$seurat_clusters,
  row.names = colnames(scRNA) # 使用colnames(scRNA)作为细胞ID
)

# 循环提取基因表达数据
for (i in 1:6) {
  gene <- gene_name[i]

  # 直接从RNA层提取表达数据（Seurat 5.0+）
  expr_data <- GetAssayData(scRNA, assay = "RNA", layer = "data")[gene, ]

  # 转换为数据框并添加到结果中
  tmp <- data.frame(gene = as.numeric(expr_data), row.names = colnames(scRNA))
  colnames(tmp) <- gene

  # 合并数据
  data_tmp <- cbind(data_tmp, tmp) # 将新列添加到data_tmp的右侧
}

# 筛选感兴趣的聚类
data_tmp <- subset(data_tmp, group %in% c(0, 1, 2, 3))

# 设置聚类水平顺序
data_tmp$group <- factor(data_tmp$group, levels = c(0, 1, 2, 3))

# 重塑数据用于绘制小提琴图
new_vio <- reshape2::melt(
  data = data_tmp,
  id.vars = "group", # 指定分组变量
  variable.name = "gene", # 基因名列名
  value.name = "exp" # 表达量列名
)

# 添加颜色分组（例如，根据聚类分组）
new_vio$fill <- ifelse(new_vio$group %in% c(0, 1), "A", "B")

# 修正后的ggplot代码
ggplot(new_vio, aes(
  x = group,
  y = exp
)) +
  geom_violin(aes(fill = fill),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#7CCD7C", "#36648B")) +
  theme_bw() +
  facet_grid(gene ~ .) + # 修正：使用实际列名"gene"而非"variable"
  xlab("") +
  ylab("") +
  theme(
    panel.grid = element_blank(),
    strip.background.x = element_blank(),
    panel.border = element_rect(linewidth = 0.5), # 修正：使用linewidth而非size
    axis.line = element_line(linewidth = 0.5), # 修正：使用linewidth而非size
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 10)
  )
ggsave(here("3_outputs", "E-MTAB-6149_outputs", "Fig2f_多个小提琴图.pdf"), width = 7, height = 10)


# 使用一个自定义的绘制多个小提琴图的R脚本stackvlion ------------------------------------------
source(here("2_src", "编程猫_src", "stackvlion.R"))
# 1. 绘制图形并赋值给对象（假设函数返回 ggplot）
p <- StackedVlnPlot(
  scRNA,
  features = gene_name,
  idents = c(0, 1, 2, 3, 4),
  cols = c("#7CCD7C", "#7CCD7C", "#36648B", "#36648B", "#36648B")
)

# 2. 保存图形（修正多余逗号，显式传图形对象）
ggsave(
  filename = here("3_outputs", "E-MTAB-6149_outputs", "Fig2f_多个小提琴图_第2种方法.pdf"),
  plot = p, # 显式指定要保存的图形对象
  height = 5,
  width = 4,
  units = "in" # 可选：明确单位（默认英寸，也可设"cm"等）
)
