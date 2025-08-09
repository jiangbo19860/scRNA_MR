
# “基因表达数据→机器学习建模→生物标志物筛选”
# 基于基因表达数据和显著基因列表，构建并评估多种机器学习模型（随机森林、SVM、XGBoost、GLM均是用于二分类任务），以实现样本组别的分类Type(Control or Tumor)?，同时分析关键基因的重要性。代码主要分为数据预处理、模型训练、模型评估和结果输出四个环节。
# 2个输入文件：合并的标准化基因表达数据1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt和显著意义的基因列表（3_outputs/汇总有显著意义的MR结果/gene.txt - geneRT）。
# 三、输出文件
# 4个评估图表
# ML_residual.pdf：模型残差的累积分布图，反映不同模型的预测误差分布。
# ML_boxplot.pdf：模型残差的箱线图，直观比较不同模型的误差离散程度。
# ML_ROC.pdf：ROC 曲线及 AUC 值，比较 4 种模型的分类性能（曲线越靠近左上角，AUC 越大，性能越好）。
# ML_importance.pdf：特征重要性图，展示各模型中对分类贡献最大的基因及其重要性分值。
# 4个重要基因列表
# importanceGene.RF.txt：随机森林模型中前 5 个最重要的基因。
# importanceGene.SVM.txt：支持向量机模型中前 5 个最重要的基因。
# importanceGene.XGB.txt：XGBoost 模型中前 5 个最重要的基因。
# importanceGene.GLM.txt：广义线性模型中前 5 个最重要的基因。
rm(list = ls())  # 清理工作环境
pacman::p_load(
  caret,      # 机器学习训练与分割
  DALEX,      # 模型解释
  ggplot2,    # 绘图
  randomForest, # 随机森林
  kernlab,    # 支持向量机
  xgboost,    # XGBoost
  pROC,       # ROC曲线分析
  glmnet,     # 正则化GLM模型
  here        # 路径管理
)

# 设置随机种子确保可重复
set.seed(123)

# 定义当天日期（用于输出文件命名）
today_date <- format(Sys.Date(), "%Y%m%d")

# 设置输入文件路径
inputFile <- here("1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt")  # 合并的基因表达数据
geneFile <- here("3_outputs/20250803_sig_genes_44.txt")  # 显著基因列表

# 读取基因表达数据（行：基因，列：样本）
data <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 筛选出列名中不包含"NA"的列（保留有效样本列）
data <- data[, !grepl("NA", colnames(data)), drop = FALSE]
cat("清理后原始数据的样本列数：", ncol(data), "\n")

# 读取基因列表文件，并提取所需的基因表达数据
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data_gene <- data[as.vector(geneRT[, 1]), , drop = FALSE]  # 提取目标基因

# 删除因基因名不匹配导致的全NA基因行（无效基因）
data_gene <- data_gene[!apply(is.na(data_gene), 1, all), , drop = FALSE]
cat("删除全NA基因后，剩余有效基因数量：", nrow(data_gene), "\n")
if (nrow(data_gene) == 0) stop("无有效基因，请检查gene.txt中的基因名是否与原始数据匹配！")

# 处理基因名中的特殊字符
row.names(data_gene) <- gsub("-", "_", row.names(data_gene))

# 转置数据：行→样本，列→基因（便于后续按样本处理）
data_transposed <- t(data_gene)
data_transposed <- as.data.frame(data_transposed)

# 删除含缺失值的样本行（确保样本数据无NA）
na_samples <- which(apply(is.na(data_transposed), 1, any))  # 定位含缺失值的样本
if (length(na_samples) > 0) {
  na_sample_ids <- rownames(data_transposed)[na_samples]
  cat("含缺失值的样本ID：", paste(na_sample_ids, collapse = ", "), "\n")
  data_clean <- data_transposed[!rownames(data_transposed) %in% na_sample_ids, , drop = FALSE]
  cat("删除缺失值样本后，剩余有效样本数量：", nrow(data_clean), "\n")
} else {
  data_clean <- data_transposed
  cat("所有样本均无缺失值，无需删除。\n")
}

# 验证清理后的数据是否无缺失值
cat("清理后样本数据中的缺失值数量：", sum(is.na(data_clean)), "\n")  # 应输出0
if (sum(is.na(data_clean)) > 0) stop("数据仍存在缺失值，请重新检查处理步骤！")

# 为清理后的样本添加组别信息，并转换为因子
data_clean$Type <- gsub("(.*)\\_(.*)", "\\2", rownames(data_clean))  # 提取组别标签
data_clean$Type <- factor(data_clean$Type, levels = c("Control", "Tumor"))  # 转换为因子，指定水平
cat("清理后样本的组别分布：", table(data_clean$Type), "\n")  # 确认组别信息正确

# 检查并删除方差为0的特征（无区分度的基因）
feature_vars_all <- setdiff(colnames(data_clean), "Type")  # 所有特征列名
feature_vars <- apply(data_clean[, feature_vars_all, drop = FALSE], 2, var)
zero_var_features <- sum(feature_vars == 0)
cat("方差为0的特征数量：", zero_var_features, "\n")

if (zero_var_features > 0) {
  valid_features <- names(feature_vars)[feature_vars > 0]
  data_clean <- data_clean[, c(valid_features, "Type")]
  cat("删除无区分度特征后，剩余有效特征数量：", length(valid_features), "\n")
  if (length(valid_features) == 0) stop("所有特征均无区分度，无法训练模型！")
}

# 划分训练集（70%）和测试集（30%）
inTrain <- createDataPartition(y = data_clean$Type, p = 0.7, list = FALSE)
train <- data_clean[inTrain, ]  # 训练集
test <- data_clean[-inTrain, ]  # 测试集
cat("训练集样本数：", nrow(train), "，测试集样本数：", nrow(test), "\n")
cat("训练集组别分布：", table(train$Type), "\n")

# 定义留一法交叉验证策略（适配小样本）
control <- trainControl(
  method = "LOOCV",  # 留一法交叉验证
  savePredictions = TRUE,
  classProbs = TRUE,  # 保留类别概率
  summaryFunction = twoClassSummary  # 二分类专用评估函数
)

# 1. 随机森林模型训练
mod_rf <- train(Type ~ ., data = train, method = "rf", trControl = control, metric = "ROC")

# 2. 支持向量机模型训练
mod_svm <- train(Type ~ ., data = train, method = "svmRadial", trControl = control, metric = "ROC")

# 3. XGBoost模型训练（特征筛选：保留最具区分度的前5个特征）
# 计算特征与类别的ANOVA相关性（筛选重要特征）
feature_anova <- apply(train[, feature_vars_all, drop = FALSE], 2, function(x) {
  tryCatch({
    model <- aov(x ~ train$Type)
    summary(model)[[1]]$`Pr(>F)`[1]  # 提取P值
  }, error = function(e) NA)
})
feature_anova <- feature_anova[!is.na(feature_anova)]  # 移除异常值

# 保留P值最小的前5个特征（最具区分度）
top_features <- names(sort(feature_anova)[1:min(5, length(feature_anova))])
train_sub <- train[, c(top_features, "Type")]  # 子集训练集（仅含重要特征）
test_sub <- test[, c(top_features, "Type")]    # 子集测试集

# 确保Type为因子
train_sub$Type <- factor(train_sub$Type, levels = c("Control", "Tumor"))
test_sub$Type <- factor(test_sub$Type, levels = c("Control", "Tumor"))

# 训练XGBoost模型
xgb_grid <- expand.grid(
  nrounds = 150,
  eta = 0.4,
  max_depth = 3,
  gamma = 0,
  colsample_bytree = c(0.6, 0.8),
  min_child_weight = 1,
  subsample = c(0.5, 0.75, 1.0)
)
mod_xgb <- train(
  Type ~ .,
  data = train_sub,
  method = "xgbTree",
  trControl = control,
  metric = "ROC",
  tuneGrid = xgb_grid
)

# 4. 改进的GLM模型训练（使用正则化解决多重共线性问题）
# 为GLM单独筛选特征（减少特征数量，避免秩亏）
# 计算特征相关性，移除高相关特征（相关系数>0.8）
cor_matrix <- cor(train[, feature_vars_all, drop = FALSE])
high_cor <- findCorrelation(cor_matrix, cutoff = 0.8)
if (length(high_cor) > 0) {
  glm_features <- feature_vars_all[-high_cor]
} else {
  glm_features <- feature_vars_all
}

# 进一步限制GLM特征数量（样本数的1/2）
max_glm_features <- min(floor(nrow(train)/2), length(glm_features))
if (max_glm_features < length(glm_features)) {
  # 按ANOVA P值选择最相关的特征
  glm_anova <- feature_anova[names(feature_anova) %in% glm_features]
  glm_features <- names(sort(glm_anova)[1:max_glm_features])
}

cat("GLM模型使用的特征数量：", length(glm_features), "\n")
train_glm <- train[, c(glm_features, "Type")]  # GLM专用训练集

# 使用弹性网正则化GLM（glmnet）解决多重共线性
glm_grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = seq(0.001, 0.1, by = 0.001))
mod_glm <- train(
  Type ~ .,
  data = train_glm,
  method = "glmnet",  # 使用正则化GLM
  family = "binomial",
  trControl = control,
  metric = "ROC",
  tuneGrid = glm_grid
)

# 定义预测函数（用于模型解释）
p_fun <- function(object, newdata) {
  # 确保newdata包含模型所需的所有特征
  if ("glmnet" %in% class(object) && !all(glm_features %in% colnames(newdata))) {
    # 为GLM模型准备正确的特征集
    newdata <- newdata[, intersect(colnames(newdata), c(glm_features, "Type")), drop = FALSE]
  }
  predict(object, newdata = newdata, type = "prob")[, "Tumor"]
}

yTest <- ifelse(test$Type == "Control", 0, 1)  # 测试集标签转换为0/1

# 准备GLM专用测试集
test_glm <- test[, c(glm_features, "Type")]

# 模型解释与性能评估
explainer_rf <- explain(mod_rf, label = "RF", data = test, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_svm <- explain(mod_svm, label = "SVM", data = test, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_xgb <- explain(mod_xgb, label = "XGB", data = test_sub, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_glm <- explain(mod_glm, label = "GLM", data = test_glm, y = yTest, predict_function = p_fun, verbose = FALSE)

# 模型性能对象
mp_rf <- model_performance(explainer_rf)
mp_svm <- model_performance(explainer_svm)
mp_xgb <- model_performance(explainer_xgb)
mp_glm <- model_performance(explainer_glm)

# 输出1：模型残差累积分布图
pdf(file = here("3_outputs", paste0(today_date, "_ML_residual.pdf")), width = 6, height = 6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

# 输出2：模型残差箱线图
pdf(file = here("3_outputs", paste0(today_date, "_ML_boxplot.pdf")), width = 6, height = 6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()

# 输出3：ROC曲线（带AUC值）
pred_rf <- predict(mod_rf, newdata = test, type = "prob")[, "Tumor"]
pred_svm <- predict(mod_svm, newdata = test, type = "prob")[, "Tumor"]
pred_xgb <- predict(mod_xgb, newdata = test_sub, type = "prob")[, "Tumor"]
pred_glm <- predict(mod_glm, newdata = test_glm, type = "prob")[, "Tumor"]

roc_rf <- roc(yTest, pred_rf)
roc_svm <- roc(yTest, pred_svm)
roc_xgb <- roc(yTest, pred_xgb)
roc_glm <- roc(yTest, pred_glm)

pdf(file = here("3_outputs", paste0(today_date, "_ML_ROC.pdf")), width = 5, height = 5)
plot(roc_rf, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "red", lwd = 2)
plot(roc_svm, print.auc = FALSE, legacy.axes = TRUE, col = "blue", add = TRUE, lwd = 2)
plot(roc_xgb, print.auc = FALSE, legacy.axes = TRUE, col = "green", add = TRUE, lwd = 2)
plot(roc_glm, print.auc = FALSE, legacy.axes = TRUE, col = "purple", add = TRUE, lwd = 2)
legend("bottomright",
       legend = c(
         paste0("RF: ", sprintf("%.3f", roc_rf$auc)),
         paste0("SVM: ", sprintf("%.3f", roc_svm$auc)),
         paste0("XGB: ", sprintf("%.3f", roc_xgb$auc)),
         paste0("GLM: ", sprintf("%.3f", roc_glm$auc))
       ),
       col = c("red", "blue", "green", "purple"),
       lwd = 2, bty = "n"
)
dev.off()

# 提取特征重要性（排除基线行）
extract_importance <- function(explainer, top_n = 10) {
  imp <- variable_importance(explainer, loss_function = loss_root_mean_square)
  imp <- imp[imp$variable != "_baseline_", ]
  # 按重要性排序并确保有值
  if (nrow(imp) > 0) {
    imp <- imp[order(-abs(imp$dropout_loss)), ]
    return(head(imp, top_n))
  } else {
    # 如果没有有效重要性，返回空数据框
    return(data.frame(variable = character(), dropout_loss = numeric(), label = character()))
  }
}

top_n <- 10
importance_rf <- extract_importance(explainer_rf, top_n)
importance_svm <- extract_importance(explainer_svm, top_n)
importance_xgb <- extract_importance(explainer_xgb, top_n)
importance_glm <- extract_importance(explainer_glm, top_n)

# 输出4：特征重要性图（处理可能的空数据）
pdf(file = here("3_outputs", paste0(today_date, "_ML_importance.pdf")), width = 8, height = 10)
# 检查是否所有重要性数据框都为空
if (all(sapply(list(importance_rf, importance_svm, importance_xgb, importance_glm), nrow) == 0)) {
  # 如果都为空，绘制提示信息
  plot(1:10, 1:10, type = "n", xlab = "", ylab = "", main = "No valid feature importance data")
  text(5, 5, "No valid feature importance could be calculated", cex = 1.2)
} else {
  # 只绘制有数据的模型
  plot_list <- list()
  if (nrow(importance_rf) > 0) plot_list[["RF"]] <- importance_rf
  if (nrow(importance_svm) > 0) plot_list[["SVM"]] <- importance_svm
  if (nrow(importance_xgb) > 0) plot_list[["XGB"]] <- importance_xgb
  if (nrow(importance_glm) > 0) plot_list[["GLM"]] <- importance_glm

  # 合并数据并绘图
  if (length(plot_list) > 0) {
    combined_imp <- do.call(rbind, lapply(names(plot_list), function(name) {
      df <- plot_list[[name]]
      df$model <- name
      df
    }))

    # 自定义重要性图
    ggplot(combined_imp, aes(x = reorder(variable, dropout_loss), y = dropout_loss, fill = model)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_wrap(~model, ncol = 1) +
      labs(x = "Feature", y = "Importance (Dropout Loss)", fill = "Model") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}
dev.off()

# 输出5：前5个重要基因列表
geneNum <- 5
save_importance <- function(importance_df, model_name) {
  if (nrow(importance_df) >= geneNum) {
    top_genes <- head(importance_df, geneNum)
    write.table(
      top_genes,
      file = here("3_outputs", paste0(today_date, "_importanceGene.", model_name, ".txt")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  } else {
    cat(paste0("Warning: Not enough features for ", model_name, " to save top ", geneNum, " genes\n"))
  }
}

save_importance(importance_rf, "RF")
save_importance(importance_svm, "SVM")
save_importance(importance_xgb, "XGB")
save_importance(importance_glm, "GLM")

cat("所有分析完成！结果已保存至3_outputs文件夹，文件名含日期：", today_date, "\n")
