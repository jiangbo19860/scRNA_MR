# 检查并安装必要的包
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("here")

# 清除工作空间
rm(list = ls())

# 加载所需的R包
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(here)

# 定义分析函数
analyze_subset <- function(data_subset, subset_name) {

  # 提取分组信息
  Type <- gsub("(.*)\\_(.*)", "\\2", colnames(data_subset))
  cat("=== 样本分组信息 (", subset_name, ") ===\n", sep = "")
  print(table(Type))  # 打印分组统计，确认分组是否正确（避免分组不平衡）

  # 过滤无变异的基因（所有样本表达值相同的基因，无法进行差异分析）
  # 计算每行（基因）的标准差，过滤掉标准差为0的基因
  gene_sd <- apply(data_subset, 1, sd, na.rm = TRUE)  # 计算每行标准差
  data_subset <- data_subset[gene_sd > 0, , drop = FALSE]  # 保留有变异的基因
  cat("过滤后保留的基因数量：", nrow(data_subset), "\n")

  if (nrow(data_subset) == 0) {
    cat(subset_name, "中无有效基因（所有基因无表达变异），停止分析\n\n")
    return(NULL)
  }

  # 差异分析
  sigVec <- c()
  sigGeneVec <- c()
  na_genes <- c()  # 记录无法计算p值的基因

  for (i in row.names(data_subset)) {
    # 提取该基因在两组中的表达值
    gene_data <- data_subset[i, ]
    group1 <- gene_data[Type == "Control"]
    group2 <- gene_data[Type == "Tumor"]

    # 检查两组是否都有足够的数据（至少1个样本，且避免全部相同）
    if (length(unique(group1)) <= 1 || length(unique(group2)) <= 1) {
      na_genes <- c(na_genes, i)
      next  # 跳过无法检验的基因
    }

    # 进行Wilcoxon检验
    test <- tryCatch({
      wilcox.test(gene_data ~ Type)
    }, error = function(e) {
      return(NULL)  # 捕获检验错误
    })

    # 处理检验结果
    if (is.null(test) || is.na(test$p.value)) {
      na_genes <- c(na_genes, i)
      next
    }

    pvalue <- test$p.value

    # 定义显著性符号
    Sig <- ifelse(pvalue < 0.001, "***",
                  ifelse(pvalue < 0.01, "**",
                         ifelse(pvalue < 0.05, "*", "")))

    if (pvalue < 0.05) {
      sigVec <- c(sigVec, paste0(i, Sig))
      sigGeneVec <- c(sigGeneVec, i)
    }
  }

  # 打印无法计算p值的基因信息
  if (length(na_genes) > 0) {
    cat(subset_name, "中无法计算p值的基因数量：", length(na_genes), "\n")
    # 如需查看具体基因，可取消下面一行注释
    # print(na_genes[1:5])  # 显示前5个
  }

  # 检查显著性基因
  if (length(sigGeneVec) == 0) {
    cat(subset_name, "中未发现显著性基因 (p < 0.05)\n\n")
    return(NULL)
  }

  # 保存显著性基因数据
  data_sig <- data_subset[sigGeneVec, ]
  outTab <- rbind(ID = colnames(data_sig), data_sig)
  output_file <- here("3_outputs", paste0("diffGeneExp_", subset_name, ".csv"))
  write.table(outTab, file = output_file, sep = ",", quote = FALSE, col.names = FALSE)
  cat("显著性基因数据已保存到", output_file, "\n")

  # 绘制热图
  row.names(data_sig) <- sigVec
  Type_df <- data.frame(Type = Type, row.names = colnames(data_sig))
  heatmap_file <- here("3_outputs", paste0("heatmap_", subset_name, ".pdf"))
  pdf(file = heatmap_file, width = 7, height = 4.5)
  pheatmap(data_sig,
           annotation = Type_df,
           color = colorRampPalette(c(rep("blue", 2), "white", rep("red", 2)))(100),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           scale = "row",
           show_colnames = FALSE,
           show_rownames = TRUE,
           fontsize = 7)
  dev.off()
  cat("热图已保存到", heatmap_file, "\n")

  # 绘制箱线图
  exp <- as.data.frame(t(data_sig))
  exp <- cbind(exp, Type = Type)
  data_long <- melt(exp, id.vars = c("Type"))
  colnames(data_long) <- c("Type", "Gene", "Expression")

  p <- ggboxplot(data_long, x = "Gene", y = "Expression", color = "Type",
                 xlab = "", ylab = "Gene expression",
                 legend.title = "Type", palette = c("blue", "red"),
                 add = "point", width = 0.8)
  p <- p + rotate_x_text(60)
  p1 <- p + stat_compare_means(aes(group = Type),
                               method = "wilcox.test",
                               symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                  symbols = c("***", "**", "*", " ")),
                               label = "p.signif")

  boxplot_file <- here("3_outputs", paste0("boxplot_", subset_name, ".pdf"))
  pdf(file = boxplot_file, width = 9, height = 5)
  print(p1)
  dev.off()
  cat("箱线图已保存到", boxplot_file, "\n")

  cat(subset_name, "中显著性基因数量：", length(sigGeneVec), "\n\n")
  return(p1)
}

# 读取并预处理数据
rt <- read.table(here("3_outputs/snpExp.txt"), header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]

# 转换为数值矩阵
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

# 合并重复基因
has_duplicate <- any(duplicated(rownames(data)))
if (has_duplicate) {
  cat("发现重复的基因名，将进行合并\n")
  data <- avereps(data)
}

# 分割样本
total_samples <- ncol(data)
cat("总样本数量：", total_samples, "\n")
data_first20 <- data[, 1:20, drop = FALSE]  # 前20个样本
data_last12 <- data[, (total_samples - 11):total_samples, drop = FALSE]  # 后12个样本

# 分析前20个样本
cat("===== 开始分析前20个样本 =====\n")
p1_first20 <- analyze_subset(data_first20, "first20")

# 分析后12个样本
cat("===== 开始分析后12个样本 =====\n")
p1_last12 <- analyze_subset(data_last12, "last12")

# 显示图形
if (!is.null(p1_first20)) print(p1_first20)
if (!is.null(p1_last12)) print(p1_last12)
