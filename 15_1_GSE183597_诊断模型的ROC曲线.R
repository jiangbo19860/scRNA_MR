# 肿瘤诊断预测模型的构建与评估流程
rm(list = ls())  # 清空工作空间

# 安装并加载必要的包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  caret, randomForest, kernlab, xgboost, pROC, glmnet,
  ggplot2, stringr, dplyr, sva, DMwR
)

# 设置随机种子（保证结果可重复）
set.seed(123)

# 获取当天日期
today_date <- format(Sys.Date(), "%Y%m%d")
cat("Analysis Date: ", today_date, "\n")

# 定义诊断模型标识
diagnostic_tag <- "Diagnostic_Model"

# ----------------------
# 输入文件设置
# ----------------------
gene_files <- c(
  GLM = "/Users/lijiangbo/scRNA_MR/3_outputs/20250811_importanceGene.GLM.txt",
  RF = "/Users/lijiangbo/scRNA_MR/3_outputs/20250811_importanceGene.RF.txt",
  SVM = "/Users/lijiangbo/scRNA_MR/3_outputs/20250811_importanceGene.SVM.txt",
  XGB = "/Users/lijiangbo/scRNA_MR/3_outputs/20250811_importanceGene.XGB.txt",
  Stacking = "/Users/lijiangbo/scRNA_MR/3_outputs/20250811_importanceGene.Stacking.txt"
)

expr_file <- "/Users/lijiangbo/scRNA_MR/1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt"
output_dir <- "/Users/lijiangbo/scRNA_MR/3_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------
# 自定义函数：重定向XGBoost日志输出
# ----------------------
# 替换原有的 suppress_xgb_warnings 函数
suppress_xgb_warnings <- function(expr) {
  # 使用 capture.output 捕获所有输出（包括标准输出和错误）
  capture_output <- function(expr) {
    output <- capture.output({
      result <- tryCatch({
        expr
      }, warning = function(w) {
        # 忽略警告
        NULL
      }, error = function(e) {
        # 传递错误
        stop(e)
      })
    })
    result
  }
  capture_output(expr)
}

# ----------------------
# 读取表达矩阵并优化预处理
# ----------------------
expr_data <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)
rownames(expr_data) <- gsub("-|\\.", "_", rownames(expr_data))
cat("Expression matrix dimensions: ", nrow(expr_data), " genes × ", ncol(expr_data), " samples\n", sep = "")

# 从样本名提取分组和批次信息
sample_groups <- factor(gsub(".*_", "", colnames(expr_data)), levels = c("Control", "Tumor"))
batch <- factor(gsub("(.*?)_.*", "\\1", colnames(expr_data)))  # 提取批次信息

# 1. 批次效应校正
cat("Correcting batch effects...\n")
mod <- model.matrix(~factor(sample_groups))
combat_expr <- ComBat(dat = expr_data, batch = batch, mod = mod)
expr_data <- combat_expr
# 在批次效应校正后添加以下代码

# 提取在单个批次中表达值全为0的基因（ComBat会记录这些基因名）
uniform_genes <- attr(combat_expr, "batch")$genes.uniform

# 验证是否获取到20个基因
cat("检测到", length(uniform_genes), "个在单个批次中表达值全为0的基因\n")

# 生成带日期的文件名
today_date <- format(Sys.Date(), "%Y%m%d")
output_file <- file.path(output_dir, paste0(today_date, "_uniform_expression_genes.csv"))

# 提取这些基因对应的行（原始表达矩阵中的信息）
uniform_genes_data <- expr_data[uniform_genes, , drop = FALSE]

# 保存为CSV文件（包含基因名和所有样本的表达值）
write.csv(
  uniform_genes_data,
  file = output_file,
  row.names = TRUE,  # 保留基因名作为行名
  quote = FALSE
)

cat("已将", length(uniform_genes), "个基因保存至:", output_file, "\n")


# 2. 异常样本检测与移除
cat("Detecting and removing outliers...\n")
pca <- prcomp(t(expr_data), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])
distances <- dist(pca_df) %>% as.matrix() %>% rowMeans()
outlier_threshold <- quantile(distances, 0.95)
normal_samples <- names(distances)[distances < outlier_threshold]

# 更新表达矩阵和样本分组
expr_data <- expr_data[, normal_samples, drop = FALSE]
sample_groups <- sample_groups[normal_samples]
cat("After outlier removal - Samples: ", length(normal_samples),
    " (Controls: ", sum(sample_groups == "Control"),
    ", Tumors: ", sum(sample_groups == "Tumor"), ")\n", sep = "")

# 划分训练集（70%）和测试集（30%）- 全局统一划分
inTrain <- createDataPartition(y = sample_groups, p = 0.7, list = FALSE)
train_samples <- colnames(expr_data)[inTrain]
test_samples <- colnames(expr_data)[-inTrain]
cat("Training set samples: ", length(train_samples),
    "; Test set samples: ", length(test_samples), "\n", sep = "")

# ----------------------
# 存储所有模型和结果
# ----------------------
models <- list()  # 存储基础模型
roc_results <- list()  # 存储ROC结果
youden_thresholds <- list()  # 存储Youden指数最优阈值
model_features <- list()  # 存储每个模型使用的特征集

# ----------------------
# 循环处理每个基础模型（GLM, RF, SVM, XGB）
# ----------------------
for (model_name in setdiff(names(gene_files), "Stacking")) {
  cat("\n===== Processing", model_name, "model =====\n")

  # 1. 读取基因文件
  gene_file <- gene_files[model_name]
  if (!file.exists(gene_file)) {
    warning("Gene file does not exist: ", gene_file, " - skipping this model!")
    next
  }
  gene_data <- read.table(gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  genes <- as.vector(gene_data[, 1])
  genes <- gsub("-|\\.", "_", genes)
  cat("Extracted", length(genes), "genes from gene file\n")

  # 2. 筛选基因
  valid_genes <- intersect(genes, rownames(expr_data))
  if (length(valid_genes) == 0) {
    warning(model_name, " has no matching genes - skipping this model!")
    next
  }
  cat("Matched", length(valid_genes), "genes in expression matrix\n")
  expr_filtered <- expr_data[valid_genes, , drop = FALSE]

  # 3. 数据转换并划分训练/测试集
  data_t <- as.data.frame(t(expr_filtered))
  data_t$Type <- factor(gsub(".*_", "", rownames(data_t)), levels = c("Control", "Tumor"))

  # 使用全局统一的训练/测试划分
  train_data <- data_t[train_samples, , drop = FALSE]
  test_data <- data_t[test_samples, , drop = FALSE]
  cat("Model-specific training samples: ", nrow(train_data),
      "; Test samples: ", nrow(test_data), "\n", sep = "")

  # 4. 处理类别不平衡（SMOTE）
  if (min(table(train_data$Type)) < 5) {  # 当少数类样本较少时
    cat("Applying SMOTE to handle class imbalance...\n")
    train_data <- SMOTE(Type ~ ., data = train_data, perc.over = 200, perc.under = 150)
  }

  # 5. 设置交叉验证参数
  cv_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )

  # 6. 训练模型（优化参数）
  if (model_name == "GLM") {
    model <- train(
      Type ~ .,
      data = train_data,
      method = "glm",
      family = "binomial",
      trControl = cv_control,
      metric = "ROC"
    )

  } else if (model_name == "RF") {
    # 优化RF参数
    model <- train(
      Type ~ .,
      data = train_data,
      method = "rf",
      ntree = 1000,  # 增加树的数量
      trControl = cv_control,
      metric = "ROC",
      tuneGrid = expand.grid(mtry = c(floor(sqrt(ncol(train_data)-1)) - 1,
                                      floor(sqrt(ncol(train_data)-1)),
                                      floor(sqrt(ncol(train_data)-1)) + 1))
    )

  } else if (model_name == "SVM") {
    # 优化SVM参数
    model <- train(
      Type ~ .,
      data = train_data,
      method = "svmRadial",
      trControl = cv_control,
      metric = "ROC",
      prob.model = TRUE,  # 开启概率模型
      tuneGrid = expand.grid(
        sigma = seq(0.01, 0.1, length.out = 5),
        C = seq(0.1, 10, length.out = 5)
      )
    )

  } else if (model_name == "XGB") {
    # 优化XGB参数（使用重定向抑制警告）
    model <- suppress_xgb_warnings({
      train(
        Type ~ .,
        data = train_data,
        method = "xgbDART",
        trControl = cv_control,
        metric = "ROC",
        tuneGrid = expand.grid(
          nrounds = c(100, 200, 300),
          max_depth = c(2, 3, 4),
          eta = c(0.01, 0.05, 0.1),
          gamma = 0,
          subsample = c(0.7, 0.8, 0.9),
          colsample_bytree = c(0.7, 0.8, 0.9),
          rate_drop = 0.1,
          skip_drop = 0.5,
          min_child_weight = 1
        )
      )
    })
  }

  # 7. 保存模型和特征集
  models[[model_name]] <- model
  model_features[[model_name]] <- valid_genes  # 记录当前模型使用的特征

  # 8. 在测试集上预测并评估（XGB预测同样抑制警告）
  if (model_name == "XGB") {
    test_pred <- suppress_xgb_warnings({
      predict(model, newdata = test_data, type = "prob")[, "Tumor"]
    })
  } else {
    test_pred <- predict(model, newdata = test_data, type = "prob")[, "Tumor"]
  }

  # 处理NA值
  if (any(is.na(test_pred))) {
    warning(model_name, " model predictions contain NA values - filled with 0.5")
    test_pred[is.na(test_pred)] <- 0.5
  }

  test_actual <- ifelse(test_data$Type == "Tumor", 1, 0)
  roc_obj <- roc(test_actual, test_pred) # test_actual是真实标签（0/1）
  auc_value <- round(roc_obj$auc, 3)
  cat(model_name, "model AUC: ", auc_value, "\n", sep = "")

  # 计算Youden指数最优阈值
  youden_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
  youden_thresholds[[model_name]] <- roc_obj$thresholds[youden_idx]

  # 9. 保存ROC曲线（含诊断模型标识）
  roc_pdf <- file.path(output_dir,
                       paste0(today_date, "_", diagnostic_tag, "_ROC_", model_name, ".pdf"))
  pdf(roc_pdf, width = 5, height = 5)
  plot(roc_obj,
       main = paste0(model_name, " ", diagnostic_tag, " ROC Curve"),
       print.auc = TRUE,
       col = "red",
       lwd = 2,
       legacy.axes = TRUE,
       xlab = "1 - Specificity",
       ylab = "Sensitivity")
  dev.off()
  cat("Individual model ROC curve saved to: ", roc_pdf, "\n")

  # 10. 存储结果
  roc_results[[model_name]] <- list(roc = roc_obj, auc = auc_value)
}

# ----------------------
# Stacking集成模型实现（优化版本）
# ----------------------
cat("\n===== Processing Stacking model =====\n")

# 读取Stacking模型的基因文件（用于后续分析）
stacking_gene_file <- gene_files["Stacking"]
gene_data <- read.table(stacking_gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stacking_genes <- as.vector(gene_data[, 1])
stacking_genes <- gsub("-|\\.", "_", stacking_genes)
cat("Extracted", length(stacking_genes), "genes from Stacking gene file\n")

# 为每个基础模型准备其专属特征集的数据
# 1. GLM模型专属数据
glm_genes <- model_features[["GLM"]]
expr_glm <- expr_data[glm_genes, , drop = FALSE]
data_glm <- as.data.frame(t(expr_glm))
data_glm$Type <- factor(gsub(".*_", "", rownames(data_glm)), levels = c("Control", "Tumor"))
train_glm <- data_glm[train_samples, , drop = FALSE]
test_glm <- data_glm[test_samples, , drop = FALSE]

# 2. RF模型专属数据
rf_genes <- model_features[["RF"]]
expr_rf <- expr_data[rf_genes, , drop = FALSE]
data_rf <- as.data.frame(t(expr_rf))
data_rf$Type <- factor(gsub(".*_", "", rownames(data_rf)), levels = c("Control", "Tumor"))
train_rf <- data_rf[train_samples, , drop = FALSE]
test_rf <- data_rf[test_samples, , drop = FALSE]

# 3. SVM模型专属数据
svm_genes <- model_features[["SVM"]]
expr_svm <- expr_data[svm_genes, , drop = FALSE]
data_svm <- as.data.frame(t(expr_svm))
data_svm$Type <- factor(gsub(".*_", "", rownames(data_svm)), levels = c("Control", "Tumor"))
train_svm <- data_svm[train_samples, , drop = FALSE]
test_svm <- data_svm[test_samples, , drop = FALSE]

# 4. XGB模型专属数据
xgb_genes <- model_features[["XGB"]]
expr_xgb <- expr_data[xgb_genes, , drop = FALSE]
data_xgb <- as.data.frame(t(expr_xgb))
data_xgb$Type <- factor(gsub(".*_", "", rownames(data_xgb)), levels = c("Control", "Tumor"))
train_xgb <- data_xgb[train_samples, , drop = FALSE]
test_xgb <- data_xgb[test_samples, , drop = FALSE]

cat("Stacking model training samples: ", nrow(train_glm),
    "; Test samples: ", nrow(test_glm), "\n", sep = "")

# 提取基模型在训练集上的预测概率（XGB预测抑制警告）
train_pred_glm <- predict(models[["GLM"]], newdata = train_glm, type = "prob")[, "Tumor"]
train_pred_rf <- predict(models[["RF"]], newdata = train_rf, type = "prob")[, "Tumor"]
train_pred_svm <- predict(models[["SVM"]], newdata = train_svm, type = "prob")[, "Tumor"]
train_pred_xgb <- suppress_xgb_warnings({
  predict(models[["XGB"]], newdata = train_xgb, type = "prob")[, "Tumor"]
})

# 检查训练集预测结果NA值
if (any(is.na(train_pred_svm))) {
  warning("SVM training set predictions contain NA values - filled with 0.5")
  train_pred_svm[is.na(train_pred_svm)] <- 0.5
}

# 提取基模型在测试集上的预测概率（XGB预测抑制警告）
test_pred_glm <- predict(models[["GLM"]], newdata = test_glm, type = "prob")[, "Tumor"]
test_pred_rf <- predict(models[["RF"]], newdata = test_rf, type = "prob")[, "Tumor"]
test_pred_svm <- predict(models[["SVM"]], newdata = test_svm, type = "prob")[, "Tumor"]
test_pred_xgb <- suppress_xgb_warnings({
  predict(models[["XGB"]], newdata = test_xgb, type = "prob")[, "Tumor"]
})

# 检查测试集预测结果NA值
if (any(is.na(test_pred_svm))) {
  warning("SVM test set predictions contain NA values - filled with 0.5")
  test_pred_svm[is.na(test_pred_svm)] <- 0.5
}

# 验证预测结果行数匹配
if (length(unique(c(
  length(train_pred_glm), length(train_pred_rf),
  length(train_pred_svm), length(train_pred_xgb),
  length(train_samples)
))) != 1) {
  stop("Training set prediction row counts do not match - check base models!")
}

if (length(unique(c(
  length(test_pred_glm), length(test_pred_rf),
  length(test_pred_svm), length(test_pred_xgb),
  length(test_samples)
))) != 1) {
  stop("Test set prediction row counts do not match - check base models!")
}

# 创建Stacking训练集和测试集
train_stacking <- data.frame(
  GLM = train_pred_glm,
  RF = train_pred_rf,
  SVM = train_pred_svm,
  XGB = train_pred_xgb,
  Type = train_glm$Type  # 样本标签一致，使用任意一个
)

test_stacking <- data.frame(
  GLM = test_pred_glm,
  RF = test_pred_rf,
  SVM = test_pred_svm,
  XGB = test_pred_xgb
)

# 定义元模型交叉验证参数
cv_control_stacking <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = FALSE
)

# 训练元模型（使用XGBoost作为元模型，抑制警告）
set.seed(123)
meta_model <- suppress_xgb_warnings({
  train(
    Type ~ .,
    data = train_stacking,
    method = "xgbDART",  # 使用XGBoost作为元模型
    trControl = cv_control_stacking,
    metric = "ROC",
    tuneGrid = expand.grid(
      nrounds = 100,
      max_depth = 2,
      eta = 0.1,
      gamma = 0,
      subsample = 0.8,
      colsample_bytree = 0.8,
      rate_drop = 0.1,
      skip_drop = 0.5,
      min_child_weight = 1
    )
  )
})

# 查看元模型最优参数
cat("\nStacking meta-model optimal parameters:\n")
print(meta_model$bestTune)

# 评估Stacking模型性能（预测时抑制警告）
y_test_stacking <- ifelse(test_glm$Type == "Tumor", 1, 0)  # 标签一致
meta_pred_prob <- suppress_xgb_warnings({
  predict(meta_model, newdata = test_stacking, type = "prob")[, "Tumor"]
})

# 计算AUC
roc_stacking <- roc(y_test_stacking, meta_pred_prob)
auc_stacking <- round(roc_stacking$auc, 3)
cat("\nStacking ensemble model AUC: ", auc_stacking, "\n", sep = "")

# 计算95%置信区间
ci_stacking <- ci.auc(roc_stacking)
cat("Stacking AUC (95% CI): ", auc_stacking, " [",
    round(ci_stacking[1], 3), ",", round(ci_stacking[3], 3), "]\n", sep = "")

# 计算Youden指数最优阈值
youden_stacking <- which.max(roc_stacking$sensitivities + roc_stacking$specificities - 1)
cat("Stacking performance at optimal threshold:\n")
cat("Sensitivity =", round(roc_stacking$sensitivities[youden_stacking], 3),
    ", Specificity =", round(roc_stacking$specificities[youden_stacking], 3), "\n")

# 保存Stacking结果
roc_results[["Stacking"]] <- list(roc = roc_stacking, auc = auc_stacking)
youden_thresholds[["Stacking"]] <- roc_stacking$thresholds[youden_stacking]

# 保存Stacking的ROC曲线（含诊断模型标识）
roc_pdf <- file.path(output_dir,
                     paste0(today_date, "_", diagnostic_tag, "_ROC_Stacking.pdf"))
pdf(roc_pdf, width = 5, height = 5)
plot(roc_stacking,
     main = paste0("Stacking ", diagnostic_tag, " ROC Curve"),
     print.auc = TRUE,
     col = "black",
     lwd = 2,
     lty = 2,  # 虚线区分
     legacy.axes = TRUE,
     xlab = "1 - Specificity",
     ylab = "Sensitivity")
dev.off()
cat("Stacking model ROC curve saved to: ", roc_pdf, "\n")

# ----------------------
# 所有模型性能对比
# ----------------------
cat("\n=== All models performance summary ===\n")
performance_summary <- data.frame(
  Model = names(roc_results),
  AUC = sapply(roc_results, function(x) round(x$auc, 3)),
  Sensitivity = sapply(names(roc_results), function(model) {
    idx <- which.max(roc_results[[model]]$roc$sensitivities +
                       roc_results[[model]]$roc$specificities - 1)
    round(roc_results[[model]]$roc$sensitivities[idx], 3)
  }),
  Specificity = sapply(names(roc_results), function(model) {
    idx <- which.max(roc_results[[model]]$roc$sensitivities +
                       roc_results[[model]]$roc$specificities - 1)
    round(roc_results[[model]]$roc$specificities[idx], 3)
  })
)

print(performance_summary)

# 保存性能汇总（含诊断模型标识）
summary_path <- file.path(output_dir,
                          paste0(today_date, "_", diagnostic_tag, "_model_performance_summary.txt"))
write.table(
  performance_summary,
  file = summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("Performance summary saved to: ", summary_path, "\n")

# 绘制所有模型的ROC比较图（含诊断模型标识）
compare_pdf <- file.path(output_dir,
                         paste0(today_date, "_", diagnostic_tag, "_ROC_All_Models.pdf"))
pdf(compare_pdf, width = 7, height = 6)

# 定义颜色和线条类型
colors <- c("purple", "red", "blue", "green", "black")
line_types <- c(1, 1, 1, 1, 2)  # Stacking用虚线

# 绘制第一个模型的ROC曲线
first_model <- names(roc_results)[1]
plot(roc_results[[first_model]]$roc,
     col = colors[1],
     lwd = 2,
     lty = line_types[1],
     legacy.axes = TRUE,
     main = paste0(diagnostic_tag, " ROC Curve Comparison"),
     xlab = "1 - Specificity",
     ylab = "Sensitivity")

# 叠加其他模型的ROC曲线
if (length(roc_results) > 1) {
  for (i in 2:length(roc_results)) {
    model_name <- names(roc_results)[i]
    plot(roc_results[[model_name]]$roc,
         col = colors[i],
         lwd = 2,
         lty = line_types[i],
         add = TRUE)
  }
}

# 添加图例
legend_names <- sapply(names(roc_results), function(x) {
  paste0(x, " (AUC = ", roc_results[[x]]$auc, ")")
})
legend("bottomright",
       legend = legend_names,
       col = colors,
       lwd = 2,
       lty = line_types,
       bty = "n",
       cex = 0.9)

dev.off()
cat("All models ROC comparison curve saved to: ", compare_pdf, "\n")

# 保存Stacking模型预测结果（含诊断模型标识）
stacking_result_path <- file.path(output_dir,
                                  paste0(today_date, "_", diagnostic_tag, "_stacking_predictions.txt"))
write.table(
  data.frame(
    Sample = rownames(test_glm), # 样本名
    Actual = test_glm$Type, # 实际标签（Tumor/Control）
    Stacking_Prob = meta_pred_prob, # Stacking模型预测为Tumor的概率
    Stacking_Pred = ifelse(meta_pred_prob >= youden_thresholds[["Stacking"]], "Tumor", "Control")  # 二分类结果
  ),
  file = stacking_result_path, # 保存路径
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("Stacking model predictions saved to: ", stacking_result_path, "\n")

# 打印所有模型的Youden最优阈值
cat("\n===== 各模型的Youden最优阈值 =====", "\n")
print(youden_thresholds)

# 将Youden最优阈值转换为数据框（便于保存）
threshold_df <- data.frame(
  Model = names(youden_thresholds),
  Youden_Optimal_Threshold = unlist(youden_thresholds),
  row.names = NULL
)

# 定义保存路径（含诊断模型标识和日期）
threshold_path <- file.path(
  output_dir,
  paste0(today_date, "_", diagnostic_tag, "_youden_thresholds.txt")
)

# 保存到本地
write.table(
  threshold_df,
  file = threshold_path,
  sep = "\t",        # 制表符分隔
  quote = FALSE,     # 不添加引号
  row.names = FALSE, # 不保留行名
  col.names = TRUE   # 保留列名
)

cat("Youden最优阈值已保存至：", threshold_path, "\n")

cat("\n===== All models processing completed =====", "\n")
