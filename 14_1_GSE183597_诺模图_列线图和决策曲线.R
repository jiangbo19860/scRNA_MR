# 清空工作空间
rm(list = ls())

# 安装并加载必要的包
pacman::p_load(here, rms, logistf, caret, ggplot2, rmda)

# 获取当天日期并格式化为字符串（如20250811）
today_date <- format(Sys.Date(), "%Y%m%d")

# 检查并创建输出文件夹（3_outputs）
output_dir <- here("3_outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("已创建输出文件夹：", output_dir, "\n")
}

# 设置输入文件路径
inputFile <- here("1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt")  # 合并后的矩阵文件
geneFile <- here("3_outputs/20250811_stable_importance_genes.Stable.txt")  # 含variable列的基因文件
# 提取基因文件名的后缀（如从"20250811_importanceGene.RF.txt"中提取"RF"）
gene_file_name <- basename(geneFile)  # 获取文件名（如"20250811_importanceGene.RF.txt"）
# 用正则表达式匹配 . 和 .txt 之间的部分（即后缀）
suffix <- gsub(".*\\.(.*)\\.txt$", "\\1", gene_file_name)
# 示例："20250811_importanceGene.RF.txt" 会提取出 "RF"
cat("从基因文件名提取的后缀：", suffix, "\n")

# 读取标准化数据文件
data <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE,
                   row.names = 1, stringsAsFactors = FALSE)
row.names(data) <- gsub("-", "_", row.names(data))  # 统一基因名格式（替换横线）
# # 删除前3列（保留第4列及以后）
# data <- data[, -(1:3), drop = FALSE]  # drop=FALSE确保结果仍是数据框
#
# # 覆盖保存原文件（注意：此操作会永久替换原文件，请谨慎！）
# write.table(data,
#             file = inputFile,  # 保存路径与原文件一致，实现覆盖
#             sep = "\t",        # 保持制表符分隔
#             row.names = TRUE,  # 保留行名（基因名）
#             col.names = TRUE,  # 保留列名（样本名）
#             quote = FALSE)     # 不添加引号
#
# cat("已删除前3列并覆盖原文件：", inputFile, "\n")
# 读取重要基因文件（从variable列提取基因名）
geneRT <- read.table(geneFile, header = TRUE, sep = "\t", check.names = FALSE,
                     fill = TRUE, stringsAsFactors = FALSE)

# 检查是否存在variable列
if(!"variable" %in% colnames(geneRT)) {
  stop("基因文件中未找到名为'variable'的列，请检查文件格式！")
}

# 提取并处理variable列中的基因名
gene_names <- as.vector(geneRT$variable)
gene_names <- gsub("-", "_", gene_names)       # 统一格式（替换横线）
gene_names <- gsub("\\.", "_", gene_names)     # 统一格式（替换点号）
gene_names <- trimws(gene_names)               # 去除前后空格
cat("从variable列提取的基因名（处理后）：\n")
print(gene_names)

# 检查基因名是否为空
if(length(gene_names) == 0 || all(is.na(gene_names))) {
  stop("从variable列提取的基因名为空，请检查文件内容！")
}

# 筛选数据中存在的基因
existing_genes <- intersect(gene_names, row.names(data))
missing_genes <- setdiff(gene_names, row.names(data))
if(length(missing_genes) > 0) {
  warning(paste("以下基因在表达矩阵中未找到：", paste(missing_genes, collapse = ", ")))
}
if(length(existing_genes) == 0) {
  stop("没有找到匹配的基因，请检查基因名是否正确！")
}

# 筛选数据并转置（样本为行，基因为列）
data <- data[existing_genes, ]
data <- t(data)
group <- gsub("(.*)\\_(.*)", "\\2", row.names(data))  # 提取样本分组（假设行名格式为"样本ID_分组"）
rt <- cbind(as.data.frame(data), Type = group)  # 合并基因数据和分组信息
# 1. 查看当前Type的所有水平（检查是否有多余分组）
cat("当前Type的所有水平：", unique(rt$Type), "\n")
# 明确将Type转换为二分类因子（即使水平正确，也需显式转换）
rt$Type <- factor(rt$Type, levels = c("Control", "Tumor"))  # 显式指定水平

# 验证转换结果
cat("转换后的因变量类型：", class(rt$Type), "\n")  # 应输出"factor"
cat("转换后的因变量水平：", levels(rt$Type), "\n")  # 应输出"Control" "Tumor"
# 生成datadist对象（rms包必需）
ddist <- datadist(rt)
options(datadist = "ddist")

# 构建模型公式（动态使用存在的基因）
formula_str <- paste("Type ~", paste(existing_genes, collapse = " + "))
cat("使用的模型公式：", formula_str, "\n")

# ----------------------
# 1. 构建模型
# ----------------------
# 传统逻辑回归（用于Nomogram，rms包兼容）
lrmModel <- lrm(as.formula(formula_str), data = rt, x = TRUE, y = TRUE)
# Firth回归（用于解决收敛问题，更稳健）
firthModel <- logistf(as.formula(formula_str), data = rt)

# ----------------------
# 2. 生成Nomogram列线图（基于lrm模型）并保存到3_outputs
# ----------------------
nomo_pdf_name <- here(output_dir, paste0(today_date, "_Nomo_", suffix, ".pdf"))   # 输出路径：3_outputs/日期_Nomo.pdf
pdf(nomo_pdf_name, width = 8, height = 6)
nomo <- nomogram(lrmModel, fun = plogis,
                 fun.at = c(0.0001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
                 lp = FALSE, funlabel = "Risk of Disease")
plot(nomo)
dev.off()
cat("Nomogram已保存至：", nomo_pdf_name, "\n")

# ----------------------
# 3. 生成校准曲线（基于Firth回归）并保存到3_outputs
# ----------------------
cali_pdf_name <- here(output_dir, paste0(today_date, "_Calibration_", suffix, ".pdf"))
pdf(cali_pdf_name, width = 5.5, height = 5.5)
# 使用Firth模型的预测概率进行校准
rt$pred_prob <- predict(firthModel, type = "response")
calibrate_df <- data.frame(
  predicted = rt$pred_prob,
  actual = as.numeric(rt$Type == "Tumor")  # 假设肿瘤组标签为"Tumor"
)
# 按预测概率分箱计算实际概率
calibrate_df$pred_bin <- cut(calibrate_df$predicted, breaks = seq(0, 1, by = 0.1))
calibrate_summary <- aggregate(actual ~ pred_bin + predicted,
                               data = calibrate_df,
                               FUN = function(x) c(mean = mean(x), n = length(x)))
# 绘制校准曲线
plot(calibrate_summary$predicted, calibrate_summary$actual[, "mean"],
     xlab = "Predicted Probability", ylab = "Actual Probability",
     main = "Calibration Curve (Firth Regression)",
     pch = 16, col = "blue", xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lty = 2, col = "red")  # 理想校准线
dev.off()
cat("校准曲线已保存至：", cali_pdf_name, "\n")

# ----------------------
# 4. 决策曲线分析（DCA）：交叉验证+Firth回归，使用rmda包添加成本-收益比
# ----------------------
# 安装并加载rmda包（如需）
if (!require("rmda")) install.packages("rmda")
library(rmda)

rt$Type_bin <- ifelse(rt$Type == "Control", 0, 1)  # 转换为0/1（Control=0，病例=1）

# 设置5折交叉验证
set.seed(123)  # 保证结果可重复
cv_folds <- createFolds(rt$Type_bin, k = 5, list = TRUE, returnTrain = FALSE)

# 存储每个阈值下的净获益（用于后续合并为rmda兼容格式）
all_nb <- data.frame()

# 循环每个交叉验证折（使用Firth回归拟合模型）
for (i in 1:length(cv_folds)) {
  test_idx <- cv_folds[[i]]
  train_data <- rt[-test_idx, ]
  test_data <- rt[test_idx, ]

  # 在训练集上拟合Firth回归模型
  fold_model <- logistf(as.formula(formula_str), data = train_data)

  # 在测试集上预测概率
  test_data$pred_prob <- predict(fold_model, newdata = test_data, type = "response")

  # 计算该折的净获益
  thresholds <- seq(0, 1, by = 0.01)
  nb <- data.frame(threshold = thresholds, net_benefit = NA, fold = i)

  for (t in seq_along(thresholds)) {
    threshold <- thresholds[t]
    tp <- sum(test_data$Type_bin == 1 & test_data$pred_prob >= threshold)  # 真阳性
    fp <- sum(test_data$Type_bin == 0 & test_data$pred_prob >= threshold)  # 假阳性
    n <- nrow(test_data)  # 总样本量
    nb$net_benefit[t] <- (tp / n) - (fp / n) * (threshold / (1 - threshold))  # 净获益公式
  }

  all_nb <- rbind(all_nb, nb)
}

# 计算所有折的平均净获益
mean_nb <- aggregate(net_benefit ~ threshold, data = all_nb, FUN = mean)

# ----------------------
# 关键修改：使用rmda包处理并绘图（支持cost.benefit.axis）
# ----------------------
# 1. 构造rmda兼容的决策曲线数据（需包含预测概率和真实标签）
# 合并所有测试集的预测概率和真实标签（用于rmda的decision_curve函数）
cv_preds <- data.frame()  # 存储所有交叉验证的预测结果
for (i in 1:length(cv_folds)) {
  test_idx <- cv_folds[[i]]
  train_data <- rt[-test_idx, ]
  test_data <- rt[test_idx, ]
  fold_model <- logistf(as.formula(formula_str), data = train_data)
  test_data$pred_prob <- predict(fold_model, newdata = test_data, type = "response")
  cv_preds <- rbind(cv_preds, test_data[, c("Type_bin", "pred_prob")])
}
colnames(cv_preds) <- c("actual", "predicted")  # rmda要求的列名格式

# 2. 使用rmda的decision_curve函数生成决策曲线数据
dc <- decision_curve(
  actual ~ predicted,  # 模型公式（实际标签~预测概率）
  data = cv_preds,
  family = "binomial",  # 二分类模型
  thresholds = seq(0, 1, by = 0.01),  # 阈值范围
  confidence.intervals = FALSE  # 关闭置信区间（如需可设为0.95）
)

# 3. 绘制DCA曲线，添加成本-收益比标注
dca_pdf_name <- here(output_dir, paste0(today_date, "_DCA_", suffix, ".pdf"))
pdf(dca_pdf_name, width = 6, height = 5.5)  # 适当加宽画布容纳双轴
plot_decision_curve(
  dc,
  curve.names = "Model",  # 曲线名称
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  main = "Decision Curve Analysis",
  cost.benefit.axis = TRUE,  # 自动添加成本-收益比轴
  col = "red",
  lwd = 2,
  confidence.intervals = FALSE
)
dev.off()

cat("带成本-收益比的DCA曲线已保存至：", dca_pdf_name, "\n")

# 输出所有结果文件路径
cat("\n分析完成，生成的文件：\n")
cat(nomo_pdf_name, "\n")
cat(cali_pdf_name, "\n")
cat(dca_pdf_name, "\n")
