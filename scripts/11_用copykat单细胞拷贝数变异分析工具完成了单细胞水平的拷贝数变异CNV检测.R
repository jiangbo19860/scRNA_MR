rm(list = ls())

# devtools::install_local("copykat-master")
# install.packages("copykat-master.zip", repos = NULL, type = "source")

pacman::p_load(here, Seurat, copykat, ggplot2, dplyr)
load(here("1_data", "E-MTAB-6149", "scRNA_endo.Rds"))

# 提取原始计数矩阵（copykat需要非标准化的count数据）
expr <- as.matrix(scRNA@assays$RNA@counts)

# 使用copykat算法进行拷贝数变异(CNV)分析，最终生成的t-SNE/UMAP图中：红色 / 深色区域：标记为aneuploid的细胞，代表具有 CNV 的肿瘤细胞。蓝色 / 浅色区域：标记为diploid的细胞，代表正常二倍体细胞（如基质细胞、免疫细胞）。
copykat.test <- copykat(
  rawmat = expr, # 输入：基因×细胞的原始计数矩阵
  id.type = "S", # 样本ID类型："S"表示单细胞ID
  cell.line = "no", # 非细胞系数据（适用于原代肿瘤样本）
  ngene.chr = 3, # 每个染色体至少需要3个基因用于分析
  LOW.DR = 0.01, # 过滤低检测率基因的阈值（<1%表达的基因被过滤）
  win.size = 25, # 滑动窗口大小（25个基因一组计算CNV）
  KS.cut = 0.01, # Kolmogorov-Smirnov检验的显著性阈值
  sam.name = "st", # 输出文件名前缀
  distance = "euclidean", # 细胞聚类的距离度量方法
  n.cores = 1 # 并行计算核心数（设为1避免内存溢出）
)


# 剔除未收敛的细胞 ----------------------------------------------------------------
# copykat分析可能会产生未收敛的细胞（如噪音细胞或低质量细胞）。未收敛意味着算法无法稳定推断该细胞的拷贝数状态（如染色体片段的扩增 / 缺失），其结果可能存在较大误差。
# 假设未收敛的细胞ID为"294"、"295"（根据实际编号修改）
bad_cells <- c("294", "295")
# 剔除细胞
scRNA_filtered <- scRNA[, !colnames(scRNA) %in% bad_cells]


# 保存CNV分析结果（避免重复计算）
save(copykat.test, file = "copykat_demo.Rds")

# 提取每个细胞的CNV预测结果
pred <- copykat.test$prediction
pred <- as.data.frame(pred)

# 确保预测结果与Seurat对象中的细胞顺序一致
pred <- pred[rownames(scRNA@meta.data), ]

# 将CNV预测结果添加到Seurat对象的元数据中
# copykat.pred列包含："aneuploid"（非整倍体，即肿瘤细胞）或"diploid"（二倍体，即正常细胞）
scRNA$copykat <- pred$copykat.pred

# 可视化CNV预测结果（在t-SNE/UMAP图上按CNV状态着色）
DimPlot(scRNA, group.by = "copykat") +
  labs(title = "Tumor vs. Normal Cells")

# 统计肿瘤细胞比例
table(scRNA$copykat)


# 统计和可视化每个染色体的CNV频率 -----------------------------------------------------------
# 查看 copykat.test 的结构
str(copykat.test)

# 从 CNAmat 中提取染色体信息和细胞的 CNV 值
chrom <- copykat.test$CNAmat$chrom # 染色体编号
cnv_values <- copykat.test$CNAmat[, -(1:4)] # 去掉前4列（chrom、chrompos、abspos、基因名），保留细胞的 CNV 值

# 检查维度
cat("染色体数量：", length(chrom), "\n")
cat("细胞数量：", ncol(cnv_values), "\n")

# 对每个细胞的 CNV 值判断是否 > 0.01
cnv_status <- as.data.frame(cnv_values > 0.01)

# 将染色体信息重复 ncol(cnv_values) 次（每个细胞对应一条染色体）
chrom_expanded <- rep(chrom, ncol(cnv_status))

# 将所有细胞的 CNV 状态展开为一维向量
status_expanded <- unlist(cnv_status)

# 生成染色体 × CNV 状态的频率表
cnv_freq <- table(
  Chromosome = chrom_expanded,
  CNV_Status = status_expanded
)

# 转换为数据框
cnv_freq_df <- as.data.frame(cnv_freq)

# 绘制柱状图
library(ggplot2)
ggplot(cnv_freq_df, aes(x = Chromosome, y = Freq, fill = CNV_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("FALSE" = "gray", "TRUE" = "red"),
    labels = c("Normal", "Amplification")
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "CNV Frequency by Chromosome",
    x = "Chromosome", y = "Frequency"
  )

# 保存上面的结果为Rds文件
save(scRNA, file = "scRNA_endo_copykat.Rds")
