rm(list = ls())
pacman::p_load(Seurat, dplyr, limma, stringr, GSVA, clusterProfiler, data.table, ggplot2, cowplot, patchwork, ggrepel, hdf5r, BiocParallel)
here()
load(here("3_outputs", "E-MTAB-6149_outputs", "scRNA_endo.Rds"))
load(here("1_data", "E-MTAB-6149", "scRNA_processed.Rds"))

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



# 绘制Fig2d图 ----------------------------------------------------------------
# scRNA新建group列，
scRNA$group <- ifelse(scRNA$CellFromTumor == "1", "Tumor", "Non_malignant") # R 语言中变量名不能包含连字符（-）, 所以要用_代替。
# 创建一个列名为"sample"和"group"的数据框
group_list <- data.frame(
  sample = colnames(gsva_result),
  group = scRNA$group
)
head(group_list)

# 创建设计矩阵（design matrix），设计矩阵的核心作用是将实验中的分组 / 变量信息转化为数学形式，让统计模型（如线性模型、limma 包的差异分析模型）能够理解样本间的分组关系，
# 1. 从group_list中提取分组信息，用factor转换为因子（明确组别水平）
# 2. 用~0 + ... 去除截距项，确保每个组别单独作为一列（避免多重共线性）
design <- model.matrix(~ 0 + factor(group_list$group))

# 给列命名为组别名称（如“Control”“Treatment”），方便后续识别
colnames(design) <- levels(factor(group_list$group))

# 行名设为样本名（与GSVA结果的样本名对应，确保一一匹配）
rownames(design) <- colnames(gsva_result)
# 查看设计矩阵
head(design)

contrast.matrix <- makeContrasts(
  Tumor - Non_malignant, # 定义对比组：Tumor组与Non-malignant组的差异
  levels = design # 指定设计矩阵
)

# 使用limma包的lmFit函数对通路活性矩阵进行线性模型拟合
# 输入：
#   gsva_result：通路活性矩阵（行=通路，列=细胞，值=通路活性分数）
#   design：设计矩阵（行=细胞，列=组别，值=0/1，标记细胞所属组别）
# 输出：
#   fit：包含线性模型拟合结果的对象，存储了每个通路在不同组别的系数、标准误等信息
fit <- lmFit(gsva_result, design)

# 使用contrasts.fit函数将对比矩阵应用到线性模型结果中
# 输入：
#   fit：lmFit的输出结果（包含各组别系数）
#   contrast.matrix：对比矩阵（定义要比较的组别差异，如Tumor - Non_malignant）
# 输出：
#   fit2：更新后的模型对象，存储了对比组的差异系数（如Tumor组比Non_malignant组的活性变化量）
fit2 <- contrasts.fit(fit, contrast.matrix)

# 使用eBayes函数进行 empirical Bayes 统计检验（贝叶斯校正的t检验）。贝叶斯方法收缩样本方差（尤其是样本量较小时），提高检验效率，相比普通 t 检验，eBayes能更好地处理单细胞数据中 “样本量小但通路多” 的场景，减少假阳性结果。
# 输入：
#   fit2：contrasts.fit的输出结果（包含组间差异系数）
# 输出：
#   fit2：更新后的模型对象，增加了t统计量、P值、校正后P值（FDR）等统计量
fit2 <- eBayes(fit2)

# 使用topTable函数提取显著差异通路的结果
# 输入：
#   fit2：eBayes的输出结果（包含统计检验信息）
#   coef = 1：指定要提取的对比项（这里是contrast.matrix中的第一个对比，即Tumor - Non_malignant）
#   n = Inf：返回所有通路（而非默认的前100个）
#   adjust.method = "BH"：使用Benjamini-Hochberg方法校正P值（控制FDR）
#   sort.by = "P"：按P值从小到大排序（最显著的通路排在前面）
# 输出：
#   x：数据框，包含每个通路的差异统计量（如logFC、t值、P值、FDR等）
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

# 查看结果数据框x的前6行
# 作用：快速预览最显著的几个通路的差异统计量，判断分析是否合理
head(x)
# 将结果保存为CSV文件，方便后续查看和分析
write.csv(x, here("3_outputs", "E-MTAB-6149_outputs", "endo_gsva_limma.csv"), quote = FALSE)
# 创建一个数据框，包含通路ID和对应的t统计量
df <- data.frame(ID = rownames(x), score = x$t)
# 按照t统计量的值分组（cutoff为阈值，分为3组）
cutoff <- 2
df$group <- cut(df$score,
  breaks = c(-Inf, -cutoff, cutoff, Inf),
  labels = c(1, 2, 3)
)
# 按照t统计量排序
sortdf <- df[order(df$score), ]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
# 查看排序后的数据框
head(sortdf)
# 绘制条形图，展示每个通路的t统计量
ggplot(sortdf, aes(ID, score, fill = group)) +
  geom_bar(stat = "identity") + # 使用条形图展示t统计量
  coord_flip() + # 翻转坐标轴，使通路名称在y轴
  scale_fill_manual(values = c("palegreen3", "snow3", "dodgerblue4"), guide = "none") + # 自定义颜色
  geom_hline(yintercept = c(-cutoff, cutoff), color = "white", linetype = 2, size = 0.3) + # 添加水平线
  geom_text(data = subset(df, score < 0), aes(x = ID, y = 0.1, label = ID, color = group), size = 3, hjust = 0) + # 添加负值标签
  geom_text(data = subset(df, score > 0), aes(x = ID, y = -0.1, label = ID, color = group), size = 3, hjust = 1) + # 添加正值标签
  labs(x = "Pathway", y = "t-statistic", title = "Differential Pathway Activity in Endothelial Cells") + # 添加标题和坐标轴标签
  theme_minimal() + # 使用简洁主题
  theme(
    axis.text.x = element_text(hjust = 1), # 保留原x轴文本设置
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank() # 去除次网格线
  )



# 编程猫给出Fig2d的ggplot部分的代码 -----------------------------------------------------------
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
      y = 0.5,
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
      y = -0.5,
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
  theme(panel.border = element_rect(linewidth = 0.6)) + # 边框粗细
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) # 去除y轴


# 保存图形到文件
ggsave(here("3_outputs", "E-MTAB-6149_outputs", "Fig2d_endo_gsva_limma_plot.pdf"), width = 7, height = 10)
