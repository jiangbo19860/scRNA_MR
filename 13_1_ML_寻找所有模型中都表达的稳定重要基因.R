# 利用差异基因的表达数据构建机器学习模型，区分样本的 “Control”（对照组）和 “Tumor”（肿瘤组），评估不同模型的性能，并筛选出对分类有重要贡献的稳定基因。
# 2个输入文件：基因表达矩阵：1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt，和差异基因列表3_outputs/20250810_original_diffGeneExp_p05_26genes.txt

rm(list = ls())
# 加载所需的R包
pacman::p_load(
  caret, randomForest, kernlab, xgboost, pROC, glmnet,
  DALEX, DALEXtra, ggplot2, here, dplyr, tibble, gridExtra
)

# 设置随机种子，确保结果可重现
set.seed(123)
# 获取当前日期，用于命名输出文件
today_date <- format(Sys.Date(), "%Y%m%d")
cat("分析日期：", today_date, "\n")


### 1. 数据输入与预处理 ###
# 定义输入文件路径
inputFile <- here("1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt")
diff_gene_file <- here("3_outputs/20250810_original_diffGeneExp_p05_26genes.txt")

# 读取差异基因文件
diff_genes <- read.delim(diff_gene_file, header = TRUE, sep = "\t", check.names = FALSE)
sig_genes <- diff_genes$ID
cat("显著基因总数：", length(sig_genes), "（前6个：", paste(head(sig_genes), collapse = ", "), "）\n")

# 保存显著基因列表
geneFile <- here("3_outputs", paste0(today_date, "_extracted_sig_genes.txt"))
write.table(
  data.frame(Gene = sig_genes),
  file = geneFile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 读取基因表达数据
expr_data <- read.table(
  inputFile,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = NULL,
  stringsAsFactors = FALSE
)
cat("原始表达数据维度：", nrow(expr_data), "个基因 ×", ncol(expr_data)-1, "个样本\n")


### 2. 重复基因处理（仅合并重复行） ###
# 识别重复基因
colnames(expr_data)[1] <- "gene_name"
duplicate_genes <- unique(expr_data$gene_name[duplicated(expr_data$gene_name)])
cat("重复基因数量：", length(duplicate_genes), "\n")

# 合并重复基因（取平均值）
duplicate_rows <- expr_data[expr_data$gene_name %in% duplicate_genes, ]
merged_duplicates <- duplicate_rows %>%
  group_by(gene_name) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  as.data.frame()

# 合并去重后的基因和非重复基因
non_duplicate_rows <- expr_data[!expr_data$gene_name %in% duplicate_genes, ]
expr_unique <- bind_rows(merged_duplicates, non_duplicate_rows) %>% as.data.frame()

# 处理行名和列名
rownames(expr_unique) <- expr_unique$gene_name
expr_unique <- expr_unique[, -1, drop = FALSE]
cat("去重后表达数据维度：", nrow(expr_unique), "个基因 ×", ncol(expr_unique), "个样本\n")

# 标准化基因名称（去除特殊字符并转为大写）
rownames(expr_unique) <- rownames(expr_unique) %>%
  gsub("-", "_", .) %>% gsub("\\.", "_", .) %>% toupper(.)
sig_genes <- sig_genes %>% gsub("-", "_", .) %>% gsub("\\.", "_", .) %>% toupper(.)


### 3. 样本筛选与显著基因匹配 ###
# 筛选有效样本（去除含NA的样本）
valid_samples <- !grepl("NA", colnames(expr_unique))
expr_filtered <- expr_unique[, valid_samples, drop = FALSE]
cat("筛选后样本数：", ncol(expr_filtered), "\n")

# 找到显著基因与表达数据中的重叠基因
existing_genes <- intersect(sig_genes, rownames(expr_filtered))
cat("显著基因与表达数据的重叠数量：", length(existing_genes), "\n")

# 如果没有重叠基因则停止分析
if (length(existing_genes) == 0) {
  stop(paste0("无匹配的显著基因！\n表达数据基因示例：", paste(head(rownames(expr_filtered)), collapse = ", "), "\n显著基因示例：", paste(head(sig_genes), collapse = ", ")))
}

# 提取用于建模的显著基因表达数据
expr_sig <- expr_filtered[existing_genes, , drop = FALSE]
cat("最终用于建模的显著基因数量：", nrow(expr_sig), "\n")


### 4. 数据格式转换 ###
# 转置表达矩阵，使样本为行，基因为列
data_samples <- as.data.frame(t(expr_sig))
# 从样本名中提取组别信息（Control或Tumor）
data_samples$Type <- gsub(".+_([^_]+)$", "\\1", rownames(data_samples))
# 检查组别提取是否正常
if (!all(unique(data_samples$Type) %in% c("Control", "Tumor"))) {
  warning("部分样本组别提取异常：", paste(unique(data_samples$Type), collapse = ", "))
}
# 将组别转换为因子类型
data_samples$Type <- factor(data_samples$Type, levels = c("Control", "Tumor"))
cat("样本组别分布：", paste(names(table(data_samples$Type)), table(data_samples$Type), sep = ":", collapse = " | "), "\n")


### 5. 特征预处理（仅去除低方差特征） ###
# 定义特征列（排除组别列）
feature_cols <- setdiff(colnames(data_samples), "Type")
# 计算各特征的方差
feature_vars <- apply(data_samples[, feature_cols, drop = FALSE], 2, var, na.rm = TRUE)
# 统计方差为0的特征数量
zero_var <- sum(feature_vars == 0)

# 如果存在方差为0的特征，则移除它们
if (zero_var > 0) {
  data_samples <- data_samples[, c(names(feature_vars)[feature_vars > 0], "Type"), drop = FALSE]
  feature_cols <- setdiff(colnames(data_samples), "Type")
  cat("删除", zero_var, "个方差为0的特征，剩余特征数：", length(feature_cols), "\n")
}

# 如果所有特征都被移除，则停止分析
if (length(feature_cols) == 0) stop("所有特征方差均为0，无法建模！")

### 5. 特征预处理（修改高相关特征去除阈值） ###
# 定义特征列（排除组别列）
feature_cols <- setdiff(colnames(data_samples), "Type")

# 计算各特征的方差
feature_vars <- apply(data_samples[, feature_cols, drop = FALSE], 2, var, na.rm = TRUE)
zero_var <- sum(feature_vars == 0)
if (zero_var > 0) {
  data_samples <- data_samples[, c(names(feature_vars)[feature_vars > 0], "Type"), drop = FALSE]
  feature_cols <- setdiff(colnames(data_samples), "Type")
  cat("删除", zero_var, "个方差为0的特征，剩余特征数：", length(feature_cols), "\n")
}

# 去除高相关特征（相关系数>0.7）
if (length(feature_cols) > 1) {  # 至少2个特征才需计算相关性
  cor_mat <- cor(data_samples[, feature_cols, drop = FALSE])  # 计算相关系数矩阵
  # 找到相关系数>0.7的特征（cutoff=0.7）
  high_cor <- findCorrelation(cor_mat, cutoff = 0.7)  # 核心修改：阈值从0.8改为0.7

  if (length(high_cor) > 0) {
    # 保留非高相关特征
    keep_cols <- feature_cols[-high_cor]
    data_samples <- data_samples[, c(keep_cols, "Type"), drop = FALSE]
    feature_cols <- keep_cols
    cat("删除", length(high_cor), "个高相关特征（相关系数>0.7），剩余特征数：", length(feature_cols), "\n")
  }
}

# 如果所有特征都被移除，则停止分析（原有代码）
if (length(feature_cols) == 0) stop("所有特征方差均为0，无法建模！")

### 6. 数据集划分与交叉验证 ###
# 划分训练集（70%）和测试集（30%）
inTrain <- createDataPartition(
  y = data_samples$Type,
  p = 0.7,
  list = FALSE,
  times = 1
)
train_data <- data_samples[inTrain, , drop = FALSE]
test_data <- data_samples[-inTrain, , drop = FALSE]
cat("训练集样本数：", nrow(train_data), "；测试集样本数：", nrow(test_data), "\n")

# 定义交叉验证参数
cv_control <- trainControl(
  method = "cv",          # 5折交叉验证
  number = 5,
  savePredictions = TRUE, # 保存预测结果
  classProbs = TRUE,      # 保存类别概率
  summaryFunction = twoClassSummary, # 二分类问题的评价函数
  sampling = "up",        # 向上采样处理不平衡数据
  verboseIter = FALSE     # 不输出迭代过程
)

### 7. 模型训练（XGBoost使用全部特征） ###
### 7.1 随机森林模型训练
cat("训练随机森林模型...\n")
set.seed(123)  # 设置随机种子确保可重复性
mod_rf <- train(
  Type ~ .,               # 公式：以Type为因变量，其他为自变量
  data = train_data,      # 训练数据
  method = "rf",          # 模型类型：随机森林
  trControl = cv_control, # 交叉验证配置
  metric = "ROC",         # 优化指标：ROC曲线下面积
  ntree = 500             # 树的数量
)


### 7.2 支持向量机模型训练
cat("训练支持向量机模型...\n")
set.seed(123)
mod_svm <- train(
  Type ~ .,
  data = train_data,
  method = "svmRadial",   # 径向基核函数的支持向量机
  trControl = cv_control,
  metric = "ROC",
  tuneLength = 5          # 参数搜索长度
)


### 7.3 XGBoost模型训练（使用全部特征）
cat("训练XGBoost模型（使用全部特征）...\n")
train_xgb <- train_data  # XGBoost训练数据（全部特征）
test_xgb <- test_data    # XGBoost测试数据（与训练集特征一致）

set.seed(123)
mod_xgb <- train(
  Type ~ .,
  data = train_xgb,
  method = "xgbTree",     # XGBoost树模型
  trControl = cv_control,
  metric = "ROC",
  tuneGrid = expand.grid( # 参数网格
    nrounds = 100,        # 迭代次数
    eta = 0.3,            # 学习率
    max_depth = 3,        # 树的最大深度
    gamma = 0,            # 分裂所需的最小损失减少值
    colsample_bytree = 0.7, # 每棵树使用的特征比例
    min_child_weight = 1, # 叶节点的最小样本权重和
    subsample = 0.8       # 每棵树使用的样本比例
  )
)


### 7.4 正则化GLM（广义线性模型）训练
cat("训练GLM模型...\n")
# 计算特征相关系数矩阵
cor_mat <- cor(train_data[, feature_cols, drop = FALSE])
# 识别高相关特征（相关系数>0.8）
high_cor <- findCorrelation(cor_mat, cutoff = 0.8)
# 确定GLM使用的特征（去除高相关特征）
glm_features <- if (length(high_cor) > 0) feature_cols[-high_cor] else feature_cols
# 限制最大特征数（不超过训练样本数的1/2）
max_glm_feat <- min(floor(nrow(train_data)/2), length(glm_features))

# 如果特征数超过限制，使用ANOVA p值筛选
if (max_glm_feat < length(glm_features)) {
  # 创建ANOVA公式
  anova_formula <- as.formula(paste("Type ~", paste(glm_features, collapse = " + ")))
  # 执行ANOVA分析
  anova_model <- aov(anova_formula, data = train_data)
  # 提取p值
  anova_p <- summary(anova_model)[[1]][, "Pr(>F)"]
  # 按p值排序并选择前max_glm_feat个特征
  glm_features <- names(sort(anova_p)[1:max_glm_feat])
}

# 输出GLM使用的特征数量
cat("GLM模型使用的特征数：", length(glm_features), "\n")
# 构建GLM的训练集和测试集
train_glm <- train_data[, c(glm_features, "Type")]
test_glm <- test_data[, c(glm_features, "Type")]

set.seed(123)
mod_glm <- train(
  Type ~ .,
  data = train_glm,
  method = "glmnet",      # 弹性网正则化GLM
  family = "binomial",    # 二分类问题使用二项分布
  trControl = cv_control,
  metric = "ROC",
  tuneGrid = expand.grid( # 正则化参数网格
    alpha = c(0, 0.5, 1), # 弹性网混合系数（0=岭回归，1=LASSO）
    lambda = seq(0.001, 0.1, length.out = 10) # 正则化强度
  )
)


### 8. 模型评估 ###
# 生成测试集的实际标签（1=Tumor，0=Control）
y_test <- as.numeric(test_data$Type == "Tumor")
cat("实际标签示例：", head(y_test), "（0/1）\n")

# 定义预测概率函数（处理不同模型的特征一致性）
predict_prob <- function(model, newdata) {
  model_feats <- setdiff(colnames(model$trainingData), ".outcome")
  newdata <- newdata[, intersect(colnames(newdata), model_feats), drop = FALSE]
  as.numeric(predict(model, newdata = newdata, type = "prob")[, "Tumor"])
}

# 计算各模型在测试集上的预测概率
pred_rf <- predict_prob(mod_rf, test_data)
pred_svm <- predict_prob(mod_svm, test_data)
pred_xgb <- predict_prob(mod_xgb, test_xgb)
pred_glm <- predict_prob(mod_glm, test_glm)

# 计算ROC曲线和AUC值
roc_rf <- roc(y_test, pred_rf)
roc_svm <- roc(y_test, pred_svm)
roc_xgb <- roc(y_test, pred_xgb)
roc_glm <- roc(y_test, pred_glm)

# 输出AUC值
cat("RF的AUC：", roc_rf$auc, "\n")
cat("SVM的AUC：", roc_svm$auc, "\n")
cat("XGB的AUC：", roc_xgb$auc, "\n")
cat("GLM的AUC：", roc_glm$auc, "\n")


### 9. 模型性能深度分析 ###
# （1）计算AUC的95%置信区间
library(pROC)  # 确保pROC包已加载

ci_rf <- ci.auc(roc_rf)
ci_svm <- ci.auc(roc_svm)
ci_xgb <- ci.auc(roc_xgb)
ci_glm <- ci.auc(roc_glm)

cat("\n=== AUC置信区间分析 ===\n")
cat("RF的AUC（95%CI）：", roc_rf$auc, " [", round(ci_rf[1], 3), ",", round(ci_rf[3], 3), "]\n")
cat("SVM的AUC（95%CI）：", roc_svm$auc, " [", round(ci_svm[1], 3), ",", round(ci_svm[3], 3), "]\n")
cat("XGB的AUC（95%CI）：", roc_xgb$auc, " [", round(ci_xgb[1], 3), ",", round(ci_xgb[3], 3), "]\n")
cat("GLM的AUC（95%CI）：", roc_glm$auc, " [", round(ci_glm[1], 3), ",", round(ci_glm[3], 3), "]\n")


# （2）比较模型稳定性（交叉验证ROC标准差）
cat("\n=== 模型稳定性分析 ===\n")
cv_sd_rf <- sd(mod_rf$resample$ROC)
cv_sd_svm <- sd(mod_svm$resample$ROC)
cv_sd_xgb <- sd(mod_xgb$resample$ROC)
cv_sd_glm <- sd(mod_glm$resample$ROC)

cat("各模型交叉验证ROC的标准差（越小越稳定）：\n")
cat("RF：", round(cv_sd_rf, 4), "；SVM：", round(cv_sd_svm, 4),
    "；XGB：", round(cv_sd_xgb, 4), "；GLM：", round(cv_sd_glm, 4), "\n")


# （3）临床实用性指标（Youden指数最优阈值）
cat("\n=== 临床实用性分析 ===\n")
# 计算各模型的Youden指数最优阈值
youden_rf <- which.max(roc_rf$sensitivities + roc_rf$specificities - 1)
youden_svm <- which.max(roc_svm$sensitivities + roc_svm$specificities - 1)
youden_xgb <- which.max(roc_xgb$sensitivities + roc_xgb$specificities - 1)
youden_glm <- which.max(roc_glm$sensitivities + roc_glm$specificities - 1)

# 输出最优阈值下的灵敏度和特异度
cat("最优阈值下的性能：\n")
cat("RF：灵敏度=", round(roc_rf$sensitivities[youden_rf], 3),
    "，特异度=", round(roc_rf$specificities[youden_rf], 3), "\n")
cat("SVM：灵敏度=", round(roc_svm$sensitivities[youden_svm], 3),
    "，特异度=", round(roc_svm$specificities[youden_svm], 3), "\n")
cat("XGB：灵敏度=", round(roc_xgb$sensitivities[youden_xgb], 3),
    "，特异度=", round(roc_xgb$specificities[youden_xgb], 3), "\n")
cat("GLM：灵敏度=", round(roc_glm$sensitivities[youden_glm], 3),
    "，特异度=", round(roc_glm$specificities[youden_glm], 3), "\n")


### 10. 创建模型解释器用于高级分析 ###
cat("\n创建模型解释器用于高级分析...\n")

# 确保DALEX和DALEXtra包正确加载
if (!require("DALEX")) {
  install.packages("DALEX", dependencies = TRUE, repos = "https://cloud.r-project.org/")
  library(DALEX)
}
if (!require("DALEXtra")) {
  install.packages("DALEXtra", dependencies = TRUE, repos = "https://cloud.r-project.org/")
  library(DALEXtra)
}

# 准备解释器所需的数据
data_explain <- train_data[, feature_cols, drop = FALSE]
y_explain <- as.numeric(train_data$Type == "Tumor")

# 为caret模型创建适配的预测函数（关键步骤）
# 随机森林模型
predict_rf <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
}
explainer_rf <- DALEX::explain(
  model = mod_rf,
  data = data_explain,
  y = y_explain,
  predict_function = predict_rf,  # 手动指定预测函数
  label = "Random Forest",
  type = "classification"
)

# SVM模型
predict_svm <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
}
explainer_svm <- DALEX::explain(
  model = mod_svm,
  data = data_explain,
  y = y_explain,
  predict_function = predict_svm,
  label = "SVM",
  type = "classification"
)

# XGBoost模型
predict_xgb <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
}
explainer_xgb <- DALEX::explain(
  model = mod_xgb,
  data = data_explain,
  y = y_explain,
  predict_function = predict_xgb,
  label = "XGBoost",
  type = "classification"
)

# GLM模型
predict_glm <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
}
explainer_glm <- DALEX::explain(
  model = mod_glm,
  data = data_explain[, intersect(colnames(data_explain), glm_features), drop = FALSE],
  y = y_explain,
  predict_function = predict_glm,
  label = "GLM",
  type = "classification"
)

### 定义所有输出文件路径 ###
# 残差累积分布图
residual_pdf_path <- here("3_outputs", paste0(today_date, "_ML_residual.pdf"))
# 残差箱线图
boxplot_pdf_path <- here("3_outputs", paste0(today_date, "_ML_boxplot.pdf"))
# ROC 曲线
roc_pdf_path <- here("3_outputs", paste0(today_date, "_ML_ROC.pdf"))
# 特征重要性图
importance_pdf_path <- here("3_outputs", paste0(today_date, "_ML_importance.pdf"))
# RF 重要基因列表
rf_gene_path <- here("3_outputs", paste0(today_date, "_importanceGene.RF.txt"))
# SVM 重要基因列表
svm_gene_path <- here("3_outputs", paste0(today_date, "_importanceGene.SVM.txt"))
# XGB 重要基因列表
xgb_gene_path <- here("3_outputs", paste0(today_date, "_importanceGene.XGB.txt"))
# GLM 重要基因列表
glm_gene_path <- here("3_outputs", paste0(today_date, "_importanceGene.GLM.txt"))
# 稳定重要基因列表
stable_gene_path <- here("3_outputs", paste0(today_date, "_stable_importance_genes.txt"))
# 显著基因列表路径
sig_gene_path <- geneFile


### 11. 残差分析 ###
cat("进行残差分析...\n")

# 模型性能对象
mp_rf <- model_performance(explainer_rf)
mp_svm <- model_performance(explainer_svm)
mp_xgb <- model_performance(explainer_xgb)
mp_glm <- model_performance(explainer_glm)

# 输出1：模型残差累积分布图
pdf(file = residual_pdf_path, width = 6, height = 6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

# 输出2：模型残差箱线图
pdf(file = boxplot_pdf_path, width = 6, height = 6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()


### 12. ROC曲线 ###
# 计算预测概率
pred_rf <- predict(mod_rf, newdata = test_data, type = "prob")[, "Tumor"]
pred_svm <- predict(mod_svm, newdata = test_data, type = "prob")[, "Tumor"]
pred_xgb <- predict(mod_xgb, newdata = test_xgb, type = "prob")[, "Tumor"]
pred_glm <- predict(mod_glm, newdata = test_glm, type = "prob")[, "Tumor"]

# 计算ROC
roc_rf <- roc(y_test, pred_rf)
roc_svm <- roc(y_test, pred_svm)
roc_xgb <- roc(y_test, pred_xgb)
roc_glm <- roc(y_test, pred_glm)

# 输出3：ROC曲线
pdf(file = roc_pdf_path, width = 5, height = 5)
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

### 13. 特征重要性分析（完整流程，包含函数定义）###

# 首先定义extract_importance函数（关键：必须在调用前定义）
extract_importance <- function(explainer, top_n = 20) {
  # 计算变量重要性（保留所有迭代结果，不提前截断）
  imp <- variable_importance(explainer, loss_function = loss_root_mean_square)
  # 排除基线行（非真实特征）
  imp <- imp[imp$variable != "_baseline_", ]

  # 验证数据有效性并按重要性降序排列
  if (nrow(imp) > 0) {
    # 按dropout_loss绝对值排序（确保重要性高的在前）
    imp <- imp[order(-abs(imp$dropout_loss)), ]
    # 返回前top_n行
    return(head(imp, top_n))
  } else {
    return(data.frame(variable = character(), dropout_loss = numeric(), label = character()))
  }
}

# 1. 提取各模型的原始重要性结果（Top25）
top_n_extract <- 25  # 提取每个模型的Top25重要基因
importance_rf_raw <- extract_importance(explainer_rf, top_n_extract)  # 现在函数已定义，可正常调用
importance_svm_raw <- extract_importance(explainer_svm, top_n_extract)
importance_xgb_raw <- extract_importance(explainer_xgb, top_n_extract)
importance_glm_raw <- extract_importance(explainer_glm, top_n_extract)

# 2. 定义汇总函数
summarize_importance <- function(raw_imp, top_n = 25) {
  if (nrow(raw_imp) == 0) {
    return(data.frame(variable = character(), mean_dropout_loss = numeric(), label = character(), n_iterations = integer()))
  }
  imp_summary <- raw_imp %>%
    group_by(variable, label) %>%
    summarise(
      mean_dropout_loss = mean(dropout_loss, na.rm = TRUE),
      n_iterations = n(),
      .groups = "drop"
    ) %>%
    arrange(-mean_dropout_loss)
  return(head(imp_summary, top_n))
}

# 3. 生成汇总后的重要性对象
importance_rf <- summarize_importance(importance_rf_raw, top_n = 25)
importance_svm <- summarize_importance(importance_svm_raw, top_n = 25)
importance_xgb <- summarize_importance(importance_xgb_raw, top_n = 25)
importance_glm <- summarize_importance(importance_glm_raw, top_n = 25)

# 4. 定义保存函数并保存结果
geneNum <- 25
save_importance <- function(summarized_imp, model_name, file_path) {
  if (nrow(summarized_imp) > 0) {
    top_genes <- head(summarized_imp, geneNum) %>%
      select(variable, mean_dropout_loss, n_iterations, label)
    write.table(
      top_genes,
      file = file_path,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    cat(paste0("保存", model_name, "的", nrow(top_genes), "个重要基因：", file_path, "\n"))
  } else {
    cat(paste0("Warning: ", model_name, "没有有效特征可保存\n"))
  }
}

# 5. 保存重要基因
save_importance(importance_rf, "RF", rf_gene_path)
save_importance(importance_svm, "SVM", svm_gene_path)
save_importance(importance_xgb, "XGB", xgb_gene_path)
save_importance(importance_glm, "GLM", glm_gene_path)

### 14. 筛选稳定出现的重要基因（修正读取逻辑） ###
cat("\n开始筛选稳定出现的重要基因...\n")

# 读取各模型的重要基因列表（兼容文件为空或列名缺失的情况）
read_importance <- function(file_path, model_name) {
  if (file.exists(file_path)) {
    # 读取文件，不预设列名，避免读取错误
    df <- tryCatch({
      read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      # 若读取失败（如文件损坏），返回空数据框
      return(data.frame())
    })

    # 检查数据框是否有效且包含必要列
    required_cols <- c("variable", "mean_dropout_loss")
    if (nrow(df) > 0 && all(required_cols %in% colnames(df))) {
      df$model <- model_name  # 添加模型标签
      return(df[, c(required_cols, "model")])  # 只选择必要的列
    } else {
      # 若列名缺失或数据为空，返回空数据框
      return(data.frame(variable = character(), mean_dropout_loss = numeric(), model = character()))
    }
  } else {
    # 若文件不存在，返回空数据框
    return(data.frame(variable = character(), mean_dropout_loss = numeric(), model = character()))
  }
}

# 读取各模型的重要基因
rf_imp <- read_importance(rf_gene_path, "RF")
svm_imp <- read_importance(svm_gene_path, "SVM")
xgb_imp <- read_importance(xgb_gene_path, "XGB")
glm_imp <- read_importance(glm_gene_path, "GLM")

# 合并所有模型的重要基因
all_imp <- bind_rows(rf_imp, svm_imp, xgb_imp, glm_imp)

if (nrow(all_imp) == 0) {
  cat("警告：没有可用的重要基因数据用于筛选\n")
} else {
  # 统计每个基因的出现频率和平均重要性
  gene_freq <- all_imp %>%
    group_by(variable) %>%
    summarise(
      occurrence = n_distinct(model),  # 出现的模型数量
      mean_importance = mean(mean_dropout_loss, na.rm = TRUE),  # 跨模型平均重要性
      max_importance = max(mean_dropout_loss, na.rm = TRUE)     # 最高重要性
    ) %>%
    arrange(desc(occurrence), desc(mean_importance))

  # 筛选至少5个稳定基因（逻辑不变）
  target_count <- 5
  stable_genes <- NULL
  stable_genes <- gene_freq %>% filter(occurrence == 4)
  if (nrow(stable_genes) < target_count) {
    need <- target_count - nrow(stable_genes)
    stable_genes <- bind_rows(stable_genes, gene_freq %>% filter(occurrence == 3) %>% head(need))
  }
  if (nrow(stable_genes) < target_count) {
    need <- target_count - nrow(stable_genes)
    stable_genes <- bind_rows(stable_genes, gene_freq %>% filter(occurrence == 2) %>% head(need))
  }
  if (nrow(stable_genes) < target_count) {
    need <- target_count - nrow(stable_genes)
    stable_genes <- bind_rows(stable_genes, gene_freq %>% filter(occurrence == 1) %>% arrange(desc(mean_importance)) %>% head(need))
  }

  # 输出结果
  cat("所有模型的重要基因统计（按稳定性排序）：\n")
  print(gene_freq, n = nrow(gene_freq))

  cat("\n筛选出的稳定重要基因（共", nrow(stable_genes), "个）：\n", sep = "")
  if (nrow(stable_genes) > 0) {
    print(stable_genes, n = nrow(stable_genes))
    write.table(
      stable_genes,
      file = stable_gene_path,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    cat("稳定重要基因已保存至：", stable_gene_path, "\n")
  } else {
    cat("未找到满足条件的稳定基因，请扩大特征提取范围\n")
  }
}

### 14. 筛选稳定出现的重要基因（确保至少5个） ###
cat("\n开始筛选稳定出现的重要基因...\n")

# 读取各模型的重要基因列表（使用汇总后的去重结果）
read_importance <- function(file_path, model_name) {
  if (file.exists(file_path)) {
    df <- read.delim(file_path, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      df$model <- model_name
      return(df[, c("variable", "mean_dropout_loss", "model")])  # 使用平均重要性
    }
  }
  return(data.frame(variable = character(), mean_dropout_loss = numeric(), model = character()))
}

rf_imp <- read_importance(rf_gene_path, "RF")
svm_imp <- read_importance(svm_gene_path, "SVM")
xgb_imp <- read_importance(xgb_gene_path, "XGB")
glm_imp <- read_importance(glm_gene_path, "GLM")

all_imp <- bind_rows(rf_imp, svm_imp, xgb_imp, glm_imp)

if (nrow(all_imp) == 0) {
  cat("警告：没有可用的重要基因数据用于筛选\n")
} else {
  # 统计每个基因的出现频率和平均重要性
  gene_freq <- all_imp %>%
    group_by(variable) %>%
    summarise(
      occurrence = n_distinct(model),  # 出现的模型数量（最大值为4）
      mean_importance = mean(mean_dropout_loss, na.rm = TRUE),  # 跨模型平均重要性
      max_importance = max(mean_dropout_loss, na.rm = TRUE)     # 模型中最高重要性
    ) %>%
    arrange(desc(occurrence), desc(mean_importance))  # 优先按出现次数排序

  # 逐步放宽条件筛选，确保至少得到5个基因
  target_count <- 5
  stable_genes <- NULL
  # 1. 优先筛选在所有4个模型中出现的基因
  stable_genes <- gene_freq %>% filter(occurrence == 4)
  # 2. 若不足，补充在3个模型中出现的基因
  if (nrow(stable_genes) < target_count) {
    need <- target_count - nrow(stable_genes)
    add_genes <- gene_freq %>% filter(occurrence == 3) %>% head(need)
    stable_genes <- bind_rows(stable_genes, add_genes)
  }
  # 3. 若仍不足，补充在2个模型中出现的基因（按重要性排序）
  if (nrow(stable_genes) < target_count) {
    need <- target_count - nrow(stable_genes)
    add_genes <- gene_freq %>% filter(occurrence == 2) %>% head(need)
    stable_genes <- bind_rows(stable_genes, add_genes)
  }
  # 4. 若仍不足，补充仅在1个模型中出现但重要性最高的基因
  if (nrow(stable_genes) < target_count) {
    need <- target_count - nrow(stable_genes)
    add_genes <- gene_freq %>% filter(occurrence == 1) %>% arrange(desc(mean_importance)) %>% head(need)
    stable_genes <- bind_rows(stable_genes, add_genes)
  }

  # 输出结果
  cat("所有模型的重要基因统计（按稳定性排序）：\n")
  print(gene_freq, n = nrow(gene_freq))

  cat("\n筛选出的稳定重要基因（共", nrow(stable_genes), "个）：\n", sep = "")
  if (nrow(stable_genes) > 0) {
    print(stable_genes, n = nrow(stable_genes))

    # 保存结果
    write.table(
      stable_genes,
      file = stable_gene_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    cat("稳定重要基因已保存至：", stable_gene_path, "\n")
  } else {
    cat("未找到满足条件的稳定基因，请进一步扩大top_n的值\n")
  }
}
# 打印所有输出文件的完整路径
cat("\n所有分析完成！输出文件保存路径如下：\n")
cat("1. 显著基因列表：", sig_gene_path, "\n")
cat("2. 残差累积分布图：", residual_pdf_path, "\n")
cat("3. 残差箱线图：", boxplot_pdf_path, "\n")
cat("4. ROC曲线：", roc_pdf_path, "\n")
cat("5. 特征重要性图：", importance_pdf_path, "\n")
cat("6. RF重要基因列表：", rf_gene_path, "\n")
cat("7. SVM重要基因列表：", svm_gene_path, "\n")
cat("8. XGB重要基因列表：", xgb_gene_path, "\n")
cat("9. GLM重要基因列表：", glm_gene_path, "\n")
cat("10. 稳定重要基因列表：", stable_gene_path, "\n")
