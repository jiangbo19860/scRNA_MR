# 基因表达数据多模型比较与DCA分析完整脚本: 效果还不如直接用Firth回归好。
# 功能：比较标准逻辑回归、Firth回归、随机森林、XGBoost和集成模型的临床预测性能

# 清空工作空间
rm(list = ls())

# 安装并加载必要的包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here, rms, logistf, caret, ggplot2,
  randomForest, xgboost, pROC, gridExtra
)

# ----------------------
# 1. 基础设置与文件夹准备
# ----------------------
# 获取当前日期（用于输出文件命名）
today_date <- format(Sys.Date(), "%Y%m%d")

# 检查并创建输出文件夹
output_dir <- here("3_outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("已创建输出文件夹：", output_dir, "\n")
}

# ----------------------
# 2. 数据读取与预处理
# ----------------------
# 设置输入文件路径
inputFile <- here("1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt")  # 表达矩阵
geneFile <- here("3_outputs/20250811_stable_importance_genes.Stable.txt")  # 重要基因文件

# 提取基因文件名后缀（用于输出文件命名）
gene_file_name <- basename(geneFile)
suffix <- gsub(".*\\.(.*)\\.txt$", "\\1", gene_file_name)
cat("从基因文件名提取的后缀：", suffix, "\n")

# 读取标准化表达矩阵（行：基因，列：样本）
data <- read.table(
  inputFile,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)
row.names(data) <- gsub("-", "_", row.names(data))  # 统一基因名格式（替换横线）

# 读取重要基因文件（含variable列）
geneRT <- read.table(
  geneFile,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

# 检查基因文件格式（必须包含variable列）
if (!"variable" %in% colnames(geneRT)) {
  stop("基因文件中未找到'variable'列，请检查文件格式！")
}

# 提取并标准化基因名
gene_names <- as.vector(geneRT$variable)
gene_names <- gsub("-", "_", gene_names)       # 替换横线
gene_names <- gsub("\\.", "_", gene_names)     # 替换点号
gene_names <- trimws(gene_names)               # 去除空格
cat("从variable列提取的基因名（处理后）：\n")
print(gene_names)

# 检查基因名有效性
if (length(gene_names) == 0 || all(is.na(gene_names))) {
  stop("未从variable列提取到有效基因名！")
}

# 筛选表达矩阵中存在的基因
existing_genes <- intersect(gene_names, row.names(data))
missing_genes <- setdiff(gene_names, row.names(data))
if (length(missing_genes) > 0) {
  warning(paste("以下基因在表达矩阵中未找到：", paste(missing_genes, collapse = ", ")))
}
if (length(existing_genes) == 0) {
  stop("没有匹配的基因，请检查基因名是否正确！")
}

# 数据转置（样本为行，基因为列）并添加分组信息
data <- data[existing_genes, , drop = FALSE]  # 筛选重要基因
data <- t(data)  # 转置：样本→行，基因→列
group <- gsub("(.*)\\_(.*)", "\\2", row.names(data))  # 从样本名提取分组（假设格式：ID_分组）
rt <- cbind(as.data.frame(data), Type = group)  # 合并数据与分组

# 标准化分组变量（确保二分类因子）
rt$Type <- factor(rt$Type, levels = c("Control", "Tumor"))
rt$Type_bin <- ifelse(rt$Type == "Control", 0, 1)  # 转换为0/1标签（Control=0，Tumor=1）

# 验证数据转换结果
cat("当前样本分组水平：", levels(rt$Type), "\n")
cat("分组分布：Control =", sum(rt$Type_bin == 0), ", Tumor =", sum(rt$Type_bin == 1), "\n")

# 生成rms包所需的datadist对象
ddist <- datadist(rt)
options(datadist = "ddist")

# 构建模型公式（动态包含所有存在的基因）
formula_str <- paste("Type ~", paste(existing_genes, collapse = " + "))
cat("使用的模型公式：", formula_str, "\n")

# ----------------------
# 3. 定义待比较的模型列表
# ----------------------
models <- list(
  "标准逻辑回归" = function(train_data, formula) {
    glm(as.formula(formula), data = train_data, family = binomial(link = "logit"))
  },
  "Firth回归" = function(train_data, formula) {
    logistf(as.formula(formula), data = train_data)
  },
  "随机森林" = function(train_data, formula) {
    randomForest(as.formula(formula), data = train_data, ntree = 500, importance = TRUE)
  },
  "XGBoost" = function(train_data, formula) {
    # 转换数据为XGBoost格式（特征矩阵+标签）
    x <- model.matrix(as.formula(formula), data = train_data)[, -1]  # 去除截距列
    y <- as.numeric(train_data$Type == "Tumor")  # 0/1标签
    xgb.DMatrix(data = x, label = y) %>%
      xgboost(params = list(
        objective = "binary:logistic",  # 二分类任务
        eval_metric = "logloss",
        max_depth = 3,  # 简单调参
        eta = 0.1
      ), nrounds = 100)
  }
)

# 集成模型（取前4种模型的预测概率平均值）
ensemble_model <- function(preds_list) {
  rowMeans(do.call(cbind, preds_list))  # 简单平均集成
}

# ----------------------
# 4. 交叉验证与模型评估
# ----------------------
set.seed(123)  # 保证可重复性
cv_folds <- createFolds(rt$Type_bin, k = 5, list = TRUE)  # 5折交叉验证

# 存储所有模型的净获益结果（每行：模型+阈值+净获益）
all_models_nb <- data.frame()

# 循环交叉验证折
for (fold in 1:length(cv_folds)) {
  cat("正在进行第", fold, "折交叉验证...\n")
  test_idx <- cv_folds[[fold]]
  train_data <- rt[-test_idx, ]
  test_data <- rt[test_idx, ]

  # 存储当前折各模型的预测概率
  fold_preds <- list()

  # 1. 拟合并预测4种基础模型
  for (model_name in names(models)) {
    cat("  训练", model_name, "...\n")
    # 训练模型
    if (model_name == "XGBoost") {
      # XGBoost需要单独处理输入格式
      model <- models[[model_name]](train_data, formula_str)
      # 测试集特征矩阵
      test_x <- model.matrix(as.formula(formula_str), data = test_data)[, -1]
      pred_prob <- predict(model, newdata = test_x)
    } else {
      model <- models[[model_name]](train_data, formula_str)
      # 预测概率（不同模型的预测函数略有差异）
      if (model_name == "随机森林") {
        pred_prob <- predict(model, newdata = test_data, type = "prob")[, "Tumor"]
      } else if (model_name == "Firth回归") {
        pred_prob <- predict(model, newdata = test_data, type = "response")
      } else {  # 标准逻辑回归
        pred_prob <- predict(model, newdata = test_data, type = "response")
      }
    }
    fold_preds[[model_name]] <- pred_prob
  }

  # 2. 集成模型预测（平均4种模型的概率）
  cat("  生成集成模型预测...\n")
  fold_preds[["集成模型"]] <- ensemble_model(fold_preds)

  # 3. 计算当前折各模型的净获益
  thresholds <- seq(0, 1, by = 0.01)  # 阈值范围
  for (model_name in names(fold_preds)) {
    pred_prob <- fold_preds[[model_name]]
    nb <- data.frame(
      model = model_name,
      threshold = thresholds,
      net_benefit = NA,
      fold = fold
    )
    # 计算每个阈值下的净获益
    for (t in seq_along(thresholds)) {
      threshold <- thresholds[t]
      tp <- sum(test_data$Type_bin == 1 & pred_prob >= threshold)  # 真阳性
      fp <- sum(test_data$Type_bin == 0 & pred_prob >= threshold)  # 假阳性
      n <- nrow(test_data)
      nb$net_benefit[t] <- (tp / n) - (fp / n) * (threshold / (1 - threshold))  # 净获益公式
    }
    all_models_nb <- rbind(all_models_nb, nb)
  }
}

# 计算各模型在所有折的平均净获益
mean_nb <- aggregate(net_benefit ~ model + threshold, data = all_models_nb, FUN = mean)

# ----------------------
# 5. 绘制DCA曲线比较所有模型
# ----------------------
dca_compare_pdf <- here(output_dir, paste0(today_date, "_DCA_Comparison_", suffix, ".pdf"))
pdf(dca_compare_pdf, width = 8, height = 6)

# 绘制各模型的DCA曲线
colors <- c("标准逻辑回归" = "blue", "Firth回归" = "green",
            "随机森林" = "orange", "XGBoost" = "purple", "集成模型" = "red")
ltys <- rep(1, 5)

# 初始化绘图
plot(mean_nb$threshold[mean_nb$model == "标准逻辑回归"],
     mean_nb$net_benefit[mean_nb$model == "标准逻辑回归"],
     type = "l", col = colors["标准逻辑回归"], lwd = 2,
     xlab = "Threshold Probability", ylab = "Net Benefit",
     main = "DCA Comparison of Models", xlim = c(0, 1), ylim = range(mean_nb$net_benefit))

# 添加其他模型曲线
for (model_name in setdiff(names(colors), "标准逻辑回归")) {
  lines(mean_nb$threshold[mean_nb$model == model_name],
        mean_nb$net_benefit[mean_nb$model == model_name],
        col = colors[model_name], lwd = 2, lty = ltys[which(names(colors) == model_name)])
}

# 添加参考线（全部治疗/全部不治疗）
p <- mean(rt$Type_bin == 1)  # 总体患病率
all_treat <- p - (1 - p) * (mean_nb$threshold[mean_nb$model == "标准逻辑回归"] /
                              (1 - mean_nb$threshold[mean_nb$model == "标准逻辑回归"]))
lines(mean_nb$threshold[mean_nb$model == "标准逻辑回归"], all_treat, lty = 2, col = "black")
abline(h = 0, lty = 3, col = "gray")

# 添加图例
legend("topright", legend = c(names(colors), "Treat All", "Treat None"),
       col = c(colors, "black", "gray"), lwd = 2, lty = c(ltys, 2, 3), bty = "n")

dev.off()
cat("模型比较DCA曲线已保存至：", dca_compare_pdf, "\n")

# ----------------------
# 6. 生成其他评估图表
# ----------------------
# 传统逻辑回归（用于Nomogram）
lrmModel <- lrm(as.formula(formula_str), data = rt, x = TRUE, y = TRUE)
# Firth回归（用于校准曲线）
firthModel <- logistf(as.formula(formula_str), data = rt)

# 6.1 生成Nomogram列线图
nomo_pdf_name <- here(output_dir, paste0(today_date, "_Nomo_", suffix, ".pdf"))
pdf(nomo_pdf_name, width = 8, height = 6)
nomo <- nomogram(lrmModel, fun = plogis,
                 fun.at = c(0.0001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
                 lp = FALSE, funlabel = "Risk of Disease")
plot(nomo)
dev.off()
cat("Nomogram列线图已保存至：", nomo_pdf_name, "\n")

# 6.2 生成校准曲线（基于Firth回归）
cali_pdf_name <- here(output_dir, paste0(today_date, "_Calibration_", suffix, ".pdf"))
pdf(cali_pdf_name, width = 5.5, height = 5.5)
# 使用Firth模型的预测概率进行校准
rt$pred_prob <- predict(firthModel, type = "response")
calibrate_df <- data.frame(
  predicted = rt$pred_prob,
  actual = as.numeric(rt$Type == "Tumor")
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
# 7. 输出所有结果文件路径
# ----------------------
cat("\n分析完成，生成的文件：\n")
cat(nomo_pdf_name, "\n")
cat(cali_pdf_name, "\n")
cat(dca_compare_pdf, "\n")
