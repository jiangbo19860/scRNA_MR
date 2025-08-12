### 离群点检测与模型优化：模型效果还不如处理前！！！
# 三种方法筛选异常样本：1高残差样本（预测效果差的样本）；GLM 模型的高影响样本（杠杆值、Cook 距离异常，可能扭曲模型）；马氏距离异常的离群点（特征分布偏离整体的样本）。最终剔除 8 个被至少 2 种方法标记的样本（如 X56_E162.Tn9_N.CEL_Control 等）。

## 1. 数据清洗：
# 无重复基因，无需合并；
# 剔除 3 个高相关特征（相关系数 > 0.7），最终保留 23 个特征（基因），减少特征冗余和多重共线性对模型的干扰。
## 离群点识别：通过三种方法筛选异常样本：
# 高残差样本（预测效果差的样本）；
# GLM 模型的高影响样本（杠杆值、Cook 距离异常，可能扭曲模型）；
# 马氏距离异常的离群点（特征分布偏离整体的样本）。
# 最终剔除 8 个被至少 2 种方法标记的样本（如 X56_E162.Tn9_N.CEL_Control 等）。
## 优化后模型性能：
# 剔除离群点后，样本量为 236（Control:102，Tumor:134），重新训练 RF 和 XGB：
# RF 的 AUC 从 0.8226 提升至 0.8325（性能改善）；
# XGB 的 AUC 从 0.8253 降至 0.8（可能因离群点对 XGB 影响更大）。
# 结果说明离群点会干扰模型稳定性，且不同模型对异常值的敏感度不同（RF 更稳健）。
rm(list = ls())
# 加载所需的R包
pacman::p_load(
  caret, randomForest, kernlab, xgboost, pROC, glmnet,
  DALEX, DALEXtra, ggplot2, here, dplyr, tibble, gridExtra,
  ggfortify, statmod
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


### 5. 特征预处理 ###
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
  # 找到相关系数>0.7的特征
  high_cor <- findCorrelation(cor_mat, cutoff = 0.7)

  if (length(high_cor) > 0) {
    # 保留非高相关特征
    keep_cols <- feature_cols[-high_cor]
    data_samples <- data_samples[, c(keep_cols, "Type"), drop = FALSE]
    feature_cols <- keep_cols
    cat("删除", length(high_cor), "个高相关特征（相关系数>0.7），剩余特征数：", length(feature_cols), "\n")
  }
}

# 如果所有特征都被移除，则停止分析
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


### 7. 初始模型训练 ###
### 7.1 随机森林模型训练
cat("训练随机森林模型...\n")
set.seed(123)
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


### 7.3 XGBoost模型训练
cat("训练XGBoost模型...\n")
train_xgb <- train_data  # XGBoost训练数据
test_xgb <- test_data    # XGBoost测试数据

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
    gamma = 0,            # 分裂所需的最小损失减少值（正则化项，值越大正则越强）
    colsample_bytree = 0.7, # 每棵树使用的特征比例
    min_child_weight = 1, # 叶节点的最小样本权重和（间接正则化）
    subsample = 0.8
  )
)


### 7.4 正则化GLM模型训练
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


### 8. 初始模型评估 ###
# 生成测试集的实际标签（1=Tumor，0=Control）
y_test <- as.numeric(test_data$Type == "Tumor")
cat("实际标签示例：", head(y_test), "（0/1）\n")

# 定义预测概率函数
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

# 输出初始AUC值
cat("\n=== 初始模型AUC值 ===\n")
cat("RF的AUC：", roc_rf$auc, "\n")
cat("SVM的AUC：", roc_svm$auc, "\n")
cat("XGB的AUC：", roc_xgb$auc, "\n")
cat("GLM的AUC：", roc_glm$auc, "\n")


### 9. 样本影响分析与离群点检测 ###
### 9.1 计算样本残差（识别预测效果差的样本）
# 提取样本ID
sample_ids <- rownames(test_data)
train_sample_ids <- rownames(train_data)

# 定义计算残差的函数
calc_residuals <- function(actual, predicted_prob) {
  actual_num <- as.numeric(actual)
  actual_num <- ifelse(actual_num == 1, 0, 1)
  resid_vals <- actual_num - predicted_prob  # 临时修改变量名
  data.frame(
    sample_id = names(actual),
    residual = resid_vals,  # 对应修改引用
    abs_residual = abs(resid_vals)
  )
}

# 计算各模型的训练集和测试集残差
## RF模型
# 训练集实际标签（已添加名称）
train_actual_rf <- train_data$Type
names(train_actual_rf) <- rownames(train_data)
train_pred_rf <- predict(mod_rf, newdata = train_data, type = "prob")[, "Tumor"]
residuals_train_rf <- calc_residuals(train_actual_rf, train_pred_rf)

# 测试集实际标签（新增名称赋值）
test_actual_rf <- test_data$Type
names(test_actual_rf) <- rownames(test_data)  # 关键：添加测试集样本ID名称
residuals_test_rf <- calc_residuals(test_actual_rf, pred_rf)

all_residuals_rf <- rbind(
  residuals_train_rf %>% mutate(data_type = "train", model = "RF"),
  residuals_test_rf %>% mutate(data_type = "test", model = "RF")
)

## XGB模型
# 训练集实际标签（新增名称赋值）
train_actual_xgb <- train_xgb$Type
names(train_actual_xgb) <- rownames(train_xgb)  # 关键：添加XGB训练集样本ID名称
train_pred_xgb <- predict(mod_xgb, newdata = train_xgb, type = "prob")[, "Tumor"]
residuals_train_xgb <- calc_residuals(train_actual_xgb, train_pred_xgb)

# 测试集实际标签（新增名称赋值）
test_actual_xgb <- test_xgb$Type
names(test_actual_xgb) <- rownames(test_xgb)  # 关键：添加XGB测试集样本ID名称
residuals_test_xgb <- calc_residuals(test_actual_xgb, pred_xgb)

all_residuals_xgb <- rbind(
  residuals_train_xgb %>% mutate(data_type = "train", model = "XGB"),
  residuals_test_xgb %>% mutate(data_type = "test", model = "XGB")
)

# 合并所有模型的残差
all_residuals <- rbind(all_residuals_rf, all_residuals_xgb)

# 按样本汇总残差（取平均绝对值）
sample_residual_summary <- all_residuals %>%
  group_by(sample_id) %>%
  summarise(mean_abs_residual = mean(abs_residual, na.rm = TRUE)) %>%
  arrange(desc(mean_abs_residual))

# 筛选残差绝对值前5%的样本
threshold_residual <- quantile(sample_residual_summary$mean_abs_residual, 0.95, na.rm = TRUE)
high_residual_samples <- sample_residual_summary %>%
  filter(mean_abs_residual >= threshold_residual) %>%
  pull(sample_id)

cat("\n高残差样本数量：", length(high_residual_samples), "\n")
cat("高残差样本ID：", paste(head(high_residual_samples), collapse = ", "), if(length(high_residual_samples) > 5) "...\n" else "\n")


### 9.2 检测高影响样本（针对GLM模型）
# 提取GLM模型的核心信息
glm_core <- mod_glm$finalModel
X_glm <- model.matrix(Type ~ ., data = train_glm)[, -1]  # 设计矩阵
y_glm <- as.numeric(train_glm$Type == "Tumor")  # 响应变量

# 计算杠杆值和Cook距离
leverage <- hatvalues(lm(y_glm ~ X_glm))
cook <- cooks.distance(lm(y_glm ~ X_glm))

# 整合结果
influence_glm <- data.frame(
  sample_id = rownames(train_glm),
  leverage = leverage,
  cook_distance = cook
)

# 筛选高影响样本
avg_leverage <- mean(leverage)
high_influence_samples_glm <- influence_glm %>%
  filter(leverage > 2*avg_leverage | cook_distance > 4/nrow(train_glm)) %>%
  pull(sample_id)

cat("\nGLM模型高影响样本数量：", length(high_influence_samples_glm), "\n")
cat("GLM高影响样本ID：", paste(head(high_influence_samples_glm), collapse = ", "), if(length(high_influence_samples_glm) > 5) "...\n" else "\n")


### 9.3 离群点检测（马氏距离）
# 提取特征数据（排除标签列）
feature_data <- data_samples[, feature_cols, drop = FALSE]

# 计算马氏距离
mahalanobis_dist <- mahalanobis(
  x = feature_data,
  center = colMeans(feature_data),
  cov = cov(feature_data)
)

# 确定离群点阈值（卡方分布，99%置信区间）
df <- ncol(feature_data)
threshold_maha <- qchisq(p = 0.99, df = df)

# 筛选离群点样本
outlier_samples <- names(mahalanobis_dist)[mahalanobis_dist > threshold_maha]

cat("\n离群点样本数量：", length(outlier_samples), "\n")
cat("离群点样本ID：", paste(head(outlier_samples), collapse = ", "), if(length(outlier_samples) > 5) "...\n" else "\n")


### 9.4 综合筛选需剔除的样本
# 合并所有候选样本
candidate_samples <- unique(c(
  high_residual_samples,
  high_influence_samples_glm,
  outlier_samples
))

# 统计样本被识别的次数
sample_counts <- table(factor(c(
  high_residual_samples,
  high_influence_samples_glm,
  outlier_samples
), levels = candidate_samples))

# 筛选被至少2个指标识别的样本
to_remove <- names(sample_counts)[sample_counts >= 2]

cat("\n最终需剔除的样本数量：", length(to_remove), "\n")
cat("剔除样本ID：", paste(head(to_remove), collapse = ", "), if(length(to_remove) > 5) "...\n" else "\n")

# 保存需剔除的样本列表
outlier_file <- here("3_outputs", paste0(today_date, "_outlier_samples.txt"))
write.table(
  data.frame(sample_id = to_remove, count = as.numeric(sample_counts[to_remove])),
  file = outlier_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("剔除样本列表已保存至：", outlier_file, "\n")


### 10. 剔除样本后重新训练模型 ###
# 保留非剔除样本
if (length(to_remove) > 0) {
  data_filtered <- data_samples[!rownames(data_samples) %in% to_remove, ]
  cat("\n剔除样本后的数据维度：", nrow(data_filtered), "样本 ×", ncol(data_filtered), "特征\n")
  cat("剔除样本后的组别分布：", paste(names(table(data_filtered$Type)), table(data_filtered$Type), sep = ":", collapse = " | "), "\n")

  # 重新划分训练集和测试集
  inTrain_new <- createDataPartition(
    y = data_filtered$Type,
    p = 0.7,
    list = FALSE,
    times = 1
  )
  train_data_new <- data_filtered[inTrain_new, , drop = FALSE]
  test_data_new <- data_filtered[-inTrain_new, , drop = FALSE]
  cat("新训练集样本数：", nrow(train_data_new), "；新测试集样本数：", nrow(test_data_new), "\n")

  # 重新训练模型
  ## 随机森林
  cat("重新训练随机森林模型...\n")
  set.seed(123)
  mod_rf_new <- train(
    Type ~ .,
    data = train_data_new,
    method = "rf",
    trControl = cv_control,
    metric = "ROC",
    ntree = 500
  )

  ## XGBoost
  cat("重新训练XGBoost模型...\n")
  set.seed(123)
  mod_xgb_new <- train(
    Type ~ .,
    data = train_data_new,
    method = "xgbTree",
    trControl = cv_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      nrounds = 100,
      eta = 0.3,
      max_depth = 3,
      gamma = 0,
      colsample_bytree = 0.7,
      min_child_weight = 1,
      subsample = 0.8,       # 每棵树使用的样本比例
      lambda = c(0.01, 0.1, 1),     # L2正则化系数（重点调整：值越小正则越弱）
      alpha = 0                     # L1正则化系数（根据需要调整）
    )
  )

  # 评估新模型性能
  pred_rf_new <- predict_prob(mod_rf_new, test_data_new)
  pred_xgb_new <- predict_prob(mod_xgb_new, test_data_new)

  y_test_new <- as.numeric(test_data_new$Type == "Tumor")

  roc_rf_new <- roc(y_test_new, pred_rf_new)
  roc_xgb_new <- roc(y_test_new, pred_xgb_new)

  # 比较模型性能
  cat("\n=== 剔除样本前后模型AUC对比 ===\n")
  cat("RF模型：原始", roc_rf$auc, "→ 新模型", roc_rf_new$auc, "\n")
  cat("XGB模型：原始", roc_xgb$auc, "→ 新模型", roc_xgb_new$auc, "\n")

  # 保存性能对比结果
  performance_file <- here("3_outputs", paste0(today_date, "_model_performance_comparison.txt"))
  performance_df <- data.frame(
    Model = c("RF", "SVM", "XGB", "GLM", "RF_new", "XGB_new"),
    AUC = c(
      as.numeric(roc_rf$auc), as.numeric(roc_svm$auc),
      as.numeric(roc_xgb$auc), as.numeric(roc_glm$auc),
      as.numeric(roc_rf_new$auc), as.numeric(roc_xgb_new$auc)
    )
  )
  write.table(
    performance_df,
    file = performance_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("模型性能对比结果已保存至：", performance_file, "\n")

  # 绘制ROC曲线对比图
  roc_comparison_file <- here("3_outputs", paste0(today_date, "_roc_comparison.pdf"))
  pdf(roc_comparison_file, width = 6, height = 6)
  plot(roc_rf, print.auc = FALSE, legacy.axes = TRUE, main = "模型ROC曲线对比",
       col = "red", lwd = 2, lty = 1)
  plot(roc_rf_new, print.auc = FALSE, legacy.axes = TRUE, col = "red",
       add = TRUE, lwd = 2, lty = 2)

  plot(roc_xgb, print.auc = FALSE, legacy.axes = TRUE, col = "green",
       add = TRUE, lwd = 2, lty = 1)
  plot(roc_xgb_new, print.auc = FALSE, legacy.axes = TRUE, col = "green",
       add = TRUE, lwd = 2, lty = 2)

  legend("bottomright",
         legend = c(
           paste0("RF: ", sprintf("%.3f", roc_rf$auc)),
           paste0("RF_new: ", sprintf("%.3f", roc_rf_new$auc)),
           paste0("XGB: ", sprintf("%.3f", roc_xgb$auc)),
           paste0("XGB_new: ", sprintf("%.3f", roc_xgb_new$auc))
         ),
         col = c("red", "red", "green", "green"),
         lwd = 2, lty = c(1, 2, 1, 2), bty = "n"
  )
  dev.off()
  cat("ROC曲线对比图已保存至：", roc_comparison_file, "\n")

} else {
  cat("\n没有需要剔除的样本，使用原始模型结果。\n")
}


### 11. 模型性能深度分析与特征重要性（基于最终模型） ###
# 选择最终模型（如果有新模型则使用新模型，否则使用原始模型）
if (exists("mod_rf_new")) {
  final_rf <- mod_rf_new
  final_xgb <- mod_xgb_new
  final_train_data <- train_data_new
  final_test_data <- test_data_new
} else {
  final_rf <- mod_rf
  final_xgb <- mod_xgb
  final_train_data <- train_data
  final_test_data <- test_data
}

# 创建模型解释器
data_explain <- final_train_data[, feature_cols, drop = FALSE]
y_explain <- as.numeric(final_train_data$Type == "Tumor")

# 随机森林模型解释器
predict_rf_final <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
}
explainer_rf_final <- DALEX::explain(
  model = final_rf,
  data = data_explain,
  y = y_explain,
  predict_function = predict_rf_final,
  label = "Random Forest",
  type = "classification"
)

# XGBoost模型解释器
predict_xgb_final <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
}
explainer_xgb_final <- DALEX::explain(
  model = final_xgb,
  data = data_explain,
  y = y_explain,
  predict_function = predict_xgb_final,
  label = "XGBoost",
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
# XGB 重要基因列表
xgb_gene_path <- here("3_outputs", paste0(today_date, "_importanceGene.XGB.txt"))
# 稳定重要基因列表
stable_gene_path <- here("3_outputs", paste0(today_date, "_stable_importance_genes.txt"))
# 显著基因列表路径
sig_gene_path <- geneFile


### 12. 残差分析 ###
cat("\n进行残差分析...\n")

# 模型性能对象
mp_rf <- model_performance(explainer_rf_final)
mp_xgb <- model_performance(explainer_xgb_final)

# 输出1：模型残差累积分布图
pdf(file = residual_pdf_path, width = 6, height = 6)
p1 <- plot(mp_rf, mp_xgb)
print(p1)
dev.off()

# 输出2：模型残差箱线图
pdf(file = boxplot_pdf_path, width = 6, height = 6)
p2 <- plot(mp_rf, mp_xgb, geom = "boxplot")
print(p2)
dev.off()


### 13. ROC曲线 ###
# 计算预测概率
pred_rf_final <- predict(final_rf, newdata = final_test_data, type = "prob")[, "Tumor"]
pred_xgb_final <- predict(final_xgb, newdata = final_test_data, type = "prob")[, "Tumor"]

# 计算ROC
roc_rf_final <- roc(as.numeric(final_test_data$Type == "Tumor"), pred_rf_final)
roc_xgb_final <- roc(as.numeric(final_test_data$Type == "Tumor"), pred_xgb_final)

# 输出3：ROC曲线
pdf(file = roc_pdf_path, width = 5, height = 5)
plot(roc_rf_final, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "red", lwd = 2)
plot(roc_xgb_final, print.auc = FALSE, legacy.axes = TRUE, col = "green", add = TRUE, lwd = 2)
legend("bottomright",
       legend = c(
         paste0("RF: ", sprintf("%.3f", roc_rf_final$auc)),
         paste0("XGB: ", sprintf("%.3f", roc_xgb_final$auc))
       ),
       col = c("red", "green"),
       lwd = 2, bty = "n"
)
dev.off()


### 14. 特征重要性分析 ###

# 定义extract_importance函数
extract_importance <- function(explainer, top_n = 20) {
  imp <- variable_importance(explainer, loss_function = loss_root_mean_square)
  imp <- imp[imp$variable != "_baseline_", ]

  if (nrow(imp) > 0) {
    imp <- imp[order(-abs(imp$dropout_loss)), ]
    return(head(imp, top_n))
  } else {
    return(data.frame(variable = character(), dropout_loss = numeric(), label = character()))
  }
}

# 提取各模型的原始重要性结果（Top25）
top_n_extract <- 25
importance_rf_raw <- extract_importance(explainer_rf_final, top_n_extract)
importance_xgb_raw <- extract_importance(explainer_xgb_final, top_n_extract)

# 定义汇总函数
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

# 生成汇总后的重要性对象
importance_rf <- summarize_importance(importance_rf_raw, top_n = 25)
importance_xgb <- summarize_importance(importance_xgb_raw, top_n = 25)

# 定义保存函数并保存结果
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

# 保存重要基因
save_importance(importance_rf, "RF", rf_gene_path)
save_importance(importance_xgb, "XGB", xgb_gene_path)


### 15. 筛选稳定出现的重要基因 ###
cat("\n开始筛选稳定出现的重要基因...\n")

# 读取各模型的重要基因列表
read_importance <- function(file_path, model_name) {
  if (file.exists(file_path)) {
    df <- tryCatch({
      read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      return(data.frame())
    })

    required_cols <- c("variable", "mean_dropout_loss")
    if (nrow(df) > 0 && all(required_cols %in% colnames(df))) {
      df$model <- model_name
      return(df[, c(required_cols, "model")])
    } else {
      return(data.frame(variable = character(), mean_dropout_loss = numeric(), model = character()))
    }
  } else {
    return(data.frame(variable = character(), mean_dropout_loss = numeric(), model = character()))
  }
}

# 读取各模型的重要基因
rf_imp <- read_importance(rf_gene_path, "RF")
xgb_imp <- read_importance(xgb_gene_path, "XGB")

# 合并所有模型的重要基因
all_imp <- bind_rows(rf_imp, xgb_imp)

if (nrow(all_imp) == 0) {
  cat("警告：没有可用的重要基因数据用于筛选\n")
} else {
  # 统计每个基因的出现频率和平均重要性
  gene_freq <- all_imp %>%
    group_by(variable) %>%
    summarise(
      occurrence = n_distinct(model),
      mean_importance = mean(mean_dropout_loss, na.rm = TRUE),
      max_importance = max(mean_dropout_loss, na.rm = TRUE)
    ) %>%
    arrange(desc(occurrence), desc(mean_importance))

  # 筛选至少5个稳定基因
  target_count <- 5
  stable_genes <- NULL
  stable_genes <- gene_freq %>% filter(occurrence == 2)
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

# 打印所有输出文件的完整路径
cat("\n所有分析完成！输出文件保存路径如下：\n")
cat("1. 显著基因列表：", sig_gene_path, "\n")
cat("2. 剔除样本列表：", outlier_file, "\n")
cat("3. 模型性能对比：", if(exists("performance_file")) performance_file else "无", "\n")
cat("4. ROC曲线对比：", if(exists("roc_comparison_file")) roc_comparison_file else "无", "\n")
cat("5. 残差累积分布图：", residual_pdf_path, "\n")
cat("6. 残差箱线图：", boxplot_pdf_path, "\n")
cat("7. ROC曲线：", roc_pdf_path, "\n")
cat("8. RF重要基因列表：", rf_gene_path, "\n")
cat("9. XGB重要基因列表：", xgb_gene_path, "\n")
cat("10. 稳定重要基因列表：", stable_gene_path, "\n")
