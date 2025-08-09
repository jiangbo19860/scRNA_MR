# 对基因表达数据进行差异分析，并通过可视化手段（热图和箱线图）展示显著差异表达的基因在不同样本组（对照组 Control 和肿瘤组 Tumor）中的表达模式。输入是差异基因表达数据3_outputs/20250804_snpExp_37genes_32samples.txt。输出有3个，两个图和一个txt文件。
rm(list = ls())  # 清除工作空间
pacman::p_load(here, limma, pheatmap, reshape2, ggpubr)  # 加载所需包
here()

# --------------------------
# 1. 读取数据并预处理
# --------------------------
# 读取表达数据
rt <- read.table(here("3_outputs/20250804_snpExp_37genes_32samples.txt"),
                 header = TRUE, sep = "\t", check.names = FALSE)

# 数据格式转换
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]  # 第一列作为基因名（行名）
exp <- rt[, 2:ncol(rt)]  # 提取表达数据（排除基因名列）
data <- matrix(as.numeric(as.matrix(exp)),
               nrow = nrow(exp),
               dimnames = list(rownames(exp), colnames(exp)))

# 处理重复基因（取平均值）
has_duplicate <- any(duplicated(rownames(data)))
data <- avereps(data)
original_data <- data  # 保存原始数据用于对比
exp_full <- data  # 定义完整数据集用于绘制所有基因的箱线图


# --------------------------
# 2. 提取分组信息并统计样本量
# --------------------------
Type <- gsub("(.*)\\_(.*)", "\\2", colnames(data))  # 从样本名提取分组（Control/Tumor）
control_count <- sum(Type == "Control")
tumor_count <- sum(Type == "Tumor")
cat("原始数据样本量：对照组", control_count, "个，肿瘤组", tumor_count, "个\n")


# --------------------------
# 3. 检测平局（重复值）并打印
# --------------------------
tie_genes <- list()  # 存储存在平局的基因信息
for (gene in row.names(data)) {
  exp_values <- data[gene, ]
  if (any(duplicated(exp_values))) {  # 检查是否有重复值
    value_counts <- table(exp_values)
    tie_values <- names(value_counts)[value_counts >= 2]  # 提取重复值
    tie_genes[[gene]] <- tie_values
  }
}

# 打印平局结果
if (length(tie_genes) > 0) {
  cat("\n发现", length(tie_genes), "个基因存在平局（重复值）：\n\n")
  for (gene in names(tie_genes)) {
    cat("基因名：", gene, "\n")
    cat("平局值：", paste(tie_genes[[gene]], collapse = ", "), "\n\n")
  }
} else {
  cat("\n数据中未发现平局（所有表达量均唯一）。\n")
}

# --------------------------
# 4. 原始数据差异分析（Wilcoxon检验）
# --------------------------
# 存储原始分析结果
sigVec <- c()         # 带显著性标记的基因名
sigGeneVec <- c()     # p<0.05的基因名
for (i in row.names(data)) {
  test <- wilcox.test(data[i, ] ~ Type, exact = FALSE)  # 关闭精确p值计算（避免警告）
  pvalue <- test$p.value
  # 定义显著性标记
  Sig <- ifelse(pvalue < 0.001, "***",
                ifelse(pvalue < 0.01, "**",
                       ifelse(pvalue < 0.05, "*", "")))
  if (pvalue < 0.05) {
    sigVec <- c(sigVec, paste0(i, Sig))
    sigGeneVec <- c(sigGeneVec, i)
  }
}
sig_count <- length(sigGeneVec)  # 显著基因数量
cat("\n原始分析：p<0.05的显著基因数量：", sig_count, "\n")


# --------------------------
# 新增：绘制所有基因的箱线图
# --------------------------
# 转换完整数据格式（使用已定义的exp_full对象）
exp_full_t <- as.data.frame(t(exp_full))  # 转置完整数据（所有基因）
exp_full_t <- cbind(exp_full_t, Type = Type)  # 添加分组信息
data_melt_full <- melt(exp_full_t, id.vars = "Type")  # 重塑为长格式
colnames(data_melt_full) <- c("Type", "Gene", "Expression")  # 重命名列名

# 绘制所有基因的箱线图
p_full <- ggboxplot(data_melt_full, x = "Gene", y = "Expression", color = "Type",
                    xlab = "", ylab = "Gene expression (all genes)",
                    legend.title = "Type",
                    palette = c("green", "purple"),  # 与显著基因图颜色区分
                    add = "point",
                    width = 0.8)
p_full <- p_full + rotate_x_text(60)  # 旋转x轴标签避免重叠
p_full <- p_full + stat_compare_means(aes(group = Type),
                                      method = "wilcox.test",
                                      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                         symbols = c("***", "**", "*", " ")),
                                      label = "p.signif")

# 保存所有基因的箱线图
today_date <- format(Sys.Date(), "%Y%m%d")  # 统一日期变量
all_genes_boxplot <- paste0(today_date, "_boxplot_all_genes.pdf")
pdf(file = here("3_outputs", all_genes_boxplot), width = 12, height = 6)
print(p_full)
dev.off()
cat("所有基因箱线图已保存至：", all_genes_boxplot, "\n")


# --------------------------
# 5. 原始数据结果输出（文件+可视化）
# --------------------------

# 5.1 保存原始差异基因数据
if (sig_count > 0) {
  data_sig <- data[sigGeneVec, ]
  outTab <- rbind(ID = colnames(data_sig), data_sig)
  diff_filename <- paste0(today_date, "_original_diffGeneExp_p05_", sig_count, "genes.txt")
  write.table(outTab, file = here("3_outputs", diff_filename),
              sep = "\t", quote = FALSE, col.names = FALSE)
  row.names(data_sig) <- sigVec  # 更新行名为带标记的基因名
} else {
  diff_filename <- paste0(today_date, "_original_diffGeneExp_p05_0genes.txt")
  file.create(here("3_outputs", diff_filename))
  data_sig <- data
}

# 5.2 绘制原始数据热图（仅显著基因）
Type_df <- as.data.frame(Type, row.names = colnames(data_sig))
colnames(Type_df) <- "Group"  # 分组列名
heatmap_filename <- paste0(today_date, "_original_heatmap_p05_", sig_count, "genes.pdf")
pdf(file = here("3_outputs", heatmap_filename), width = 7, height = 4.5)
pheatmap(data_sig,
         annotation = Type_df,  # 样本分组注释
         color = colorRampPalette(c(rep("blue", 2), "white", rep("red", 2)))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize = 7)
dev.off()
cat("显著基因热图已保存至：", heatmap_filename, "\n")

# 5.3 绘制原始数据箱线图（仅显著基因）
if (sig_count > 0) {
  exp_t <- as.data.frame(t(data_sig))
  exp_t <- cbind(exp_t, Type = Type)
  data_melt <- melt(exp_t, id.vars = "Type")
  colnames(data_melt) <- c("Type", "Gene", "Expression")

  p <- ggboxplot(data_melt, x = "Gene", y = "Expression", color = "Type",
                 xlab = "", ylab = "Gene expression (significant genes)",
                 legend.title = "Group", palette = c("blue", "red"),
                 add = "point", width = 0.8) +
    rotate_x_text(60) +
    stat_compare_means(aes(group = Type),
                       method = "wilcox.test",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", " ")),
                       label = "p.signif")

  boxplot_filename <- paste0(today_date, "_original_boxplot_p05_", sig_count, "genes.pdf")
  pdf(file = here("3_outputs", boxplot_filename), width = 9, height = 5)
  print(p)
  dev.off()
  cat("显著基因箱线图已保存至：", boxplot_filename, "\n")
} else {
  boxplot_filename <- paste0(today_date, "_original_boxplot_p05_0genes.pdf")
  file.create(here("3_outputs", boxplot_filename))
  cat("无显著基因，已创建空箱线图文件：", boxplot_filename, "\n")
}


# --------------------------
# 6. 针对问题基因重新分析（剔除0值样本）
# --------------------------
# 6.1 定义存在0值平局的基因
problem_genes <- c("BLNK", "CC2D2B", "CYP2D6", "DKKL1", "EHMT2", "MYBPC2", "SLFN12L", "TNFRSF13C")
problem_genes <- intersect(problem_genes, rownames(data))  # 仅保留数据中存在的基因
cat("\n需要重新分析的问题基因数量：", length(problem_genes), "\n")

# 6.2 存储重新分析结果（对比原始p值和过滤后p值）
reanalysis_results <- data.frame(
  gene = character(),
  original_p = numeric(),
  filtered_p = numeric(),
  sig_original = character(),
  sig_filtered = character(),
  stringsAsFactors = FALSE
)

# 6.3 提取原始p值（用于对比）
for (gene in problem_genes) {
  test_original <- wilcox.test(original_data[gene, ] ~ Type, exact = FALSE)
  p_original <- test_original$p.value
  sig_original <- ifelse(p_original < 0.001, "***",
                         ifelse(p_original < 0.01, "**",
                                ifelse(p_original < 0.05, "*", "")))
  reanalysis_results <- rbind(reanalysis_results,
                              data.frame(gene = gene,
                                         original_p = p_original,
                                         filtered_p = NA,
                                         sig_original = sig_original,
                                         sig_filtered = ""))
}

# 6.4 剔除0值样本后重新分析
filtered_genes_data <- list()  # 存储过滤后的基因表达数据（用于绘图）
for (gene in problem_genes) {
  gene_exp <- original_data[gene, ]
  valid_samples <- names(gene_exp)[gene_exp != 0]  # 剔除表达量为0的样本
  if (length(valid_samples) < 2) {  # 样本量不足时跳过
    cat("基因", gene, "剔除0值后样本量不足（<2），无法分析\n")
    next
  }
  # 提取有效样本的表达量和分组
  filtered_data <- gene_exp[valid_samples]
  filtered_Type <- gsub("(.*)\\_(.*)", "\\2", valid_samples)
  if (length(unique(filtered_Type)) < 2) {  # 确保两组均有样本
    cat("基因", gene, "剔除0值后仅剩一组样本，无法分析\n")
    next
  }
  # 重新检验
  test_filtered <- wilcox.test(filtered_data ~ filtered_Type, exact = FALSE)
  p_filtered <- test_filtered$p.value
  sig_filtered <- ifelse(p_filtered < 0.001, "***",
                         ifelse(p_filtered < 0.01, "**",
                                ifelse(p_filtered < 0.05, "*", "")))
  # 更新结果表
  reanalysis_results$filtered_p[reanalysis_results$gene == gene] <- p_filtered
  reanalysis_results$sig_filtered[reanalysis_results$gene == gene] <- sig_filtered
  # 保存过滤后的基因数据（用于绘图）
  filtered_genes_data[[gene]] <- data.frame(
    Sample = names(filtered_data),
    Expression = as.numeric(filtered_data),
    Type = filtered_Type,
    Gene = gene,
    stringsAsFactors = FALSE
  )
  # 打印结果
  cat("基因", gene, "重新分析完成：原始p=", round(reanalysis_results$original_p[reanalysis_results$gene == gene], 4),
      "，过滤后p=", round(p_filtered, 4), "\n")
}

# 6.5 保存重新分析结果表
reanalysis_file <- paste0(today_date, "_reanalysis_without_zero_samples.txt")
write.table(reanalysis_results,
            file = here("3_outputs", reanalysis_file),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n重新分析结果已保存至：", reanalysis_file, "\n")


# --------------------------
# 7. 重新分析结果可视化（箱线图）
# --------------------------
if (length(filtered_genes_data) > 0) {
  # 合并过滤后的基因数据
  filtered_plot_data <- do.call(rbind, filtered_genes_data)
  # 绘制箱线图
  p_filtered <- ggboxplot(filtered_plot_data, x = "Gene", y = "Expression", color = "Type",
                          xlab = "", ylab = "Gene expression (after removing 0-values)",
                          legend.title = "Group", palette = c("green", "orange"),
                          add = "point", width = 0.8) +
    rotate_x_text(60) +
    stat_compare_means(aes(group = Type),
                       method = "wilcox.test",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", " ")),
                       label = "p.signif")

  filtered_boxplot_filename <- paste0(today_date, "_filtered_boxplot_problem_genes.pdf")
  pdf(file = here("3_outputs", filtered_boxplot_filename), width = 9, height = 5)
  print(p_filtered)
  dev.off()
  cat("过滤后基因箱线图已保存至：", filtered_boxplot_filename, "\n")
} else {
  cat("无有效过滤后基因数据，未生成重新分析箱线图\n")
}


# --------------------------
# 8. 分析总结
# --------------------------
cat("\n===== 分析总结 =====\n")
cat("1. 原始数据显著基因（p<0.05）：", sig_count, "个\n")
cat("2. 重新分析的问题基因：", length(problem_genes), "个\n")
cat("3. 所有结果已保存至 3_outputs 文件夹\n")
