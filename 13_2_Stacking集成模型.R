# 继续之前的分析流程，添加Stacking集成模型

### 15. Stacking集成模型实现 ###
cat("\n开始训练Stacking集成模型...\n")

# 提取基模型在训练集上的预测概率（作为元模型的输入特征）
train_pred_rf <- predict(mod_rf, newdata = train_data, type = "prob")[, "Tumor"]
train_pred_svm <- predict(mod_svm, newdata = train_data, type = "prob")[, "Tumor"]
train_pred_xgb <- predict(mod_xgb, newdata = train_xgb, type = "prob")[, "Tumor"]
train_pred_glm <- predict(mod_glm, newdata = train_glm, type = "prob")[, "Tumor"]

# 提取基模型在测试集上的预测概率
test_pred_rf <- predict(mod_rf, newdata = test_data, type = "prob")[, "Tumor"]
test_pred_svm <- predict(mod_svm, newdata = test_data, type = "prob")[, "Tumor"]
test_pred_xgb <- predict(mod_xgb, newdata = test_xgb, type = "prob")[, "Tumor"]
test_pred_glm <- predict(mod_glm, newdata = test_glm, type = "prob")[, "Tumor"]

# 创建用于Stacking的训练集（基模型预测结果作为特征，保留原始标签）
train_stacking <- data.frame(
  RF = train_pred_rf,
  SVM = train_pred_svm,
  XGB = train_pred_xgb,
  GLM = train_pred_glm,
  Type = train_data$Type  # 原始标签
)

# 创建用于Stacking的测试集（仅包含基模型预测结果）
test_stacking <- data.frame(
  RF = test_pred_rf,
  SVM = test_pred_svm,
  XGB = test_pred_xgb,
  GLM = test_pred_glm
)

# 定义元模型的交叉验证参数（与基模型保持一致）
cv_control_stacking <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  sampling = "up",
  verboseIter = FALSE
)

# 训练元模型（使用弹性网逻辑回归作为元模型）
set.seed(123)
meta_model <- train(
  Type ~ .,  # 以基模型预测结果为特征，预测原始标签
  data = train_stacking,
  method = "glmnet", # glmnet（弹性网逻辑回归）是线性模型，结合了岭回归和 LASSO 回归的优点，在进行特征选择的同时能对基模型的权重进行估计。
  trControl = cv_control_stacking,
  metric = "ROC",
  family = "binomial",
  tuneGrid = expand.grid(
    alpha = c(0, 0.5, 1),  # 0=岭回归, 1=LASSO, 0.5=弹性网
    lambda = seq(0.001, 0.1, length.out = 10)  # 正则化强度
  )
)

# 查看元模型的最优参数
cat("\nStacking元模型最优参数：\n")
print(meta_model$bestTune)


### 16. 评估Stacking集成模型性能 ###
# 生成测试集实际标签
y_test_stacking <- as.numeric(test_data$Type == "Tumor")

# 元模型预测概率
meta_pred_prob <- predict(meta_model, newdata = test_stacking, type = "prob")[, "Tumor"]

# 计算AUC
roc_stacking <- roc(y_test_stacking, meta_pred_prob)
cat("\nStacking集成模型的AUC：", roc_stacking$auc, "\n")

# 计算95%置信区间
ci_stacking <- ci.auc(roc_stacking)
cat("Stacking的AUC（95%CI）：", roc_stacking$auc, " [",
    round(ci_stacking[1], 3), ",", round(ci_stacking[3], 3), "]\n")

# 计算Youden指数最优阈值
youden_stacking <- which.max(roc_stacking$sensitivities + roc_stacking$specificities - 1)
cat("Stacking最优阈值下的性能：\n")
cat("灵敏度=", round(roc_stacking$sensitivities[youden_stacking], 3),
    "，特异度=", round(roc_stacking$specificities[youden_stacking], 3), "\n")


### 17. 所有模型性能对比 ###
cat("\n=== 所有模型性能汇总 ===\n")
performance_summary <- data.frame(
  Model = c("RF", "SVM", "XGB", "GLM", "Stacking"),
  AUC = c(
    as.numeric(roc_rf$auc),
    as.numeric(roc_svm$auc),
    as.numeric(roc_xgb$auc),
    as.numeric(roc_glm$auc),
    as.numeric(roc_stacking$auc)
  ),
  Sensitivity = c(
    round(roc_rf$sensitivities[youden_rf], 3),
    round(roc_svm$sensitivities[youden_svm], 3),
    round(roc_xgb$sensitivities[youden_xgb], 3),
    round(roc_glm$sensitivities[youden_glm], 3),
    round(roc_stacking$sensitivities[youden_stacking], 3)
  ),
  Specificity = c(
    round(roc_rf$specificities[youden_rf], 3),
    round(roc_svm$specificities[youden_svm], 3),
    round(roc_xgb$specificities[youden_xgb], 3),
    round(roc_glm$specificities[youden_glm], 3),
    round(roc_stacking$specificities[youden_stacking], 3)
  )
)

print(performance_summary)


### 18. 保存Stacking模型结果 ###
# 定义Stacking模型输出路径
stacking_result_path <- here("3_outputs", paste0(today_date, "_stacking_predictions.txt"))
stacking_roc_path <- here("3_outputs", paste0(today_date, "_stacking_ROC.pdf"))

# 保存预测结果
write.table(
  data.frame(
    Sample = rownames(test_data),
    Actual = test_data$Type,
    Stacking_Prob = meta_pred_prob,
    Stacking_Pred = ifelse(meta_pred_prob >= roc_stacking$thresholds[youden_stacking], "Tumor", "Control")
  ),
  file = stacking_result_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("\nStacking模型预测结果已保存至：", stacking_result_path, "\n")

# 打开PDF绘图设备，设置图片宽度6英寸、高度6英寸
pdf(file = stacking_roc_path, width = 6, height = 6)

# 绘制随机森林（RF）的ROC曲线
plot(roc_rf,
     print.auc = FALSE,                # 在图中显示AUC值
     print.auc.x = 0.6, print.auc.y = 0.2,  # 设置AUC值显示位置（x=0.6, y=0.2）
     legacy.axes = TRUE,              # 使用传统坐标轴标签（x轴：1-特异度；y轴：灵敏度）
     main = "ROC Curves (Including Stacking)",  # 图表标题
     col = "red",                     # 曲线颜色设为红色
     lwd = 2                          # 曲线线条宽度设为2（加粗）
     # 注释掉填充区域，如需可取消注释
     # auc.polygon = TRUE,
     # auc.polygon.col = adjustcolor("red", alpha = 0.1)
)

# 绘制支持向量机（SVM）的ROC曲线（添加到现有图中）
plot(roc_svm,
     print.auc = FALSE,                # 显示AUC值
     print.auc.x = 0.6, print.auc.y = 0.15,  # AUC显示位置（x=0.6, y=0.15，与RF错开）
     legacy.axes = TRUE,              # 传统坐标轴
     col = "blue",                    # 曲线颜色设为蓝色
     add = TRUE,                      # 不新建图，添加到上一张图中
     lwd = 2                          # 线条宽度2
     # 注释掉填充区域，如需可取消注释
     # auc.polygon = TRUE,
     # auc.polygon.col = adjustcolor("blue", alpha = 0.1)
)

# 绘制XGBoost（XGB）的ROC曲线（添加到现有图中）
plot(roc_xgb,
     print.auc = FALSE,                # 显示AUC值
     print.auc.x = 0.6, print.auc.y = 0.1,   # AUC显示位置（x=0.6, y=0.1）
     legacy.axes = TRUE,              # 传统坐标轴
     col = "green",                   # 曲线颜色设为绿色
     add = TRUE,                      # 添加到现有图中
     lwd = 2                          # 线条宽度2
     # 注释掉填充区域，如需可取消注释
     # auc.polygon = TRUE,
     # auc.polygon.col = adjustcolor("green", alpha = 0.1)
)

# 绘制广义线性模型（GLM）的ROC曲线（添加到现有图中）
plot(roc_glm,
     print.auc = FALSE,                # 显示AUC值
     print.auc.x = 0.6, print.auc.y = 0.05,  # AUC显示位置（x=0.6, y=0.05）
     legacy.axes = TRUE,              # 传统坐标轴
     col = "purple",                  # 曲线颜色设为紫色
     add = TRUE,                      # 添加到现有图中
     lwd = 2                          # 线条宽度2
     # 注释掉填充区域，如需可取消注释
     # auc.polygon = TRUE,
     # auc.polygon.col = adjustcolor("purple", alpha = 0.1)
)

# 绘制Stacking集成模型的ROC曲线（添加到现有图中）
plot(roc_stacking,
     print.auc = FALSE,                # 显示AUC值
     print.auc.x = 0.6, print.auc.y = 0,    # AUC显示位置（x=0.6, y=0，最底部）
     legacy.axes = TRUE,              # 传统坐标轴
     col = "black",                   # 曲线颜色设为黑色
     lwd = 2,                         # 线条宽度2
     lty = 2,                         # 线条类型设为虚线（区别于基础模型的实线）
     add = TRUE                       # 添加到现有图中
     # 注释掉填充区域，如需可取消注释
     # auc.polygon = TRUE,
     # auc.polygon.col = adjustcolor("black", alpha = 0.1)
)

# 拼接“模型名称 + AUC值”的图例标签
# 手动提取各模型的AUC值（也可从performance_summary中读取）
legend_labels <- c(
  sprintf("RF (AUC: %.3f)", as.numeric(roc_rf$auc)),
  sprintf("SVM (AUC: %.3f)", as.numeric(roc_svm$auc)),
  sprintf("XGB (AUC: %.3f)", as.numeric(roc_xgb$auc)),
  sprintf("GLM (AUC: %.3f)", as.numeric(roc_glm$auc)),
  sprintf("Stacking (AUC: %.3f)", as.numeric(roc_stacking$auc))
)

# 定义与曲线对应的颜色
legend_colors <- c("red", "blue", "green", "purple", "black")

# 添加图例（位于右下角），设置文字颜色与线条颜色一致
legend("bottomright",
       legend = legend_labels,        # 拼接后的图例标签
       col = legend_colors,           # 图例线条颜色（与曲线对应）
       text.col = legend_colors,      # 图例文字颜色（与曲线对应）
       lwd = 2,                       # 图例线条宽度（与曲线一致）
       lty = c(1, 1, 1, 1, 2),        # 图例线条类型（Stacking为虚线）
       bty = "n",                      # 去掉图例边框（更简洁）
       cex = 0.8                       # 图例文字大小（稍微缩小一点）
)

# 关闭PDF绘图设备，保存图片
dev.off()

# 打印提示信息，说明图片保存路径
cat("包含透明曲线下面积的ROC曲线已保存至：", stacking_roc_path, "\n")

# 19. 从Stacking模型中筛选重要特征基因 ###
cat("\n开始计算Stacking模型下的重要特征基因...\n")

### 步骤1：读取各基模型的特征重要性数据
# 定义读取函数（兼容不同模型的输出格式）
read_model_importance <- function(file_path) {
  if (file.exists(file_path)) {
    df <- read.delim(file_path, sep = "\t", stringsAsFactors = FALSE)
    # 确保列名正确（variable为基因名，mean_dropout_loss为重要性）
    if ("variable" %in% colnames(df) && "mean_dropout_loss" %in% colnames(df)) {
      return(df[, c("variable", "mean_dropout_loss")])
    }
  }
  return(data.frame(variable = character(), mean_dropout_loss = numeric()))
}

# 读取4个基模型的重要基因及其重要性
rf_importance <- read_model_importance(rf_gene_path) %>%
  mutate(model = "RF")
svm_importance <- read_model_importance(svm_gene_path) %>%
  mutate(model = "SVM")
xgb_importance <- read_model_importance(xgb_gene_path) %>%
  mutate(model = "XGB")
glm_importance <- read_model_importance(glm_gene_path) %>%
  mutate(model = "GLM")

# 合并所有基模型的重要性数据
all_base_importance <- bind_rows(
  rf_importance,
  svm_importance,
  xgb_importance,
  glm_importance
)

if (nrow(all_base_importance) == 0) {
  stop("无法获取基模型的特征重要性数据，请检查文件路径")
}


### 步骤2：获取元模型（Stacking）对基模型的依赖权重
# 提取glmnet元模型的系数（反映各基模型的贡献度）
meta_coef <- coef(meta_model$finalModel, meta_model$bestTune$lambda)
meta_weights <- data.frame(
  model = rownames(meta_coef)[-1],  # 排除截距项
  weight = as.numeric(meta_coef[-1, 1])  # 基模型权重
)

# 标准化权重（确保总和为1，便于加权计算）
meta_weights$weight <- abs(meta_weights$weight)  # 取绝对值（重要性与方向无关）
meta_weights$weight <- meta_weights$weight / sum(meta_weights$weight)
cat("\nStacking元模型对基模型的依赖权重：\n")
print(meta_weights)


### 步骤3：计算基因的综合重要性（基模型重要性 × 元模型权重）
gene_combined_importance <- all_base_importance %>%
  # 关联基模型权重
  left_join(meta_weights, by = "model") %>%
  # 计算加权重要性
  mutate(weighted_importance = mean_dropout_loss * weight) %>%
  # 按基因汇总（取平均加权重要性）
  group_by(variable) %>%
  summarise(
    total_importance = mean(weighted_importance, na.rm = TRUE)
  ) %>%
  # 按重要性降序排序
  arrange(desc(total_importance)) %>%
  # 过滤无效基因（排除空值或重要性为0的）
  filter(!is.na(variable), total_importance > 0)

if (nrow(gene_combined_importance) == 0) {
  stop("未计算出有效的基因重要性，请检查基模型输出")
}


### 步骤4：提取前5个重要基因
top5_genes <- head(gene_combined_importance, 5)
cat("\nStacking模型下的Top5重要特征基因：\n")
print(top5_genes)


### 步骤5：保存结果
stacking_top5_path <- here("3_outputs", paste0(today_date, "_importanceGene.Stacking.txt"))
write.table(
  top5_genes,
  file = stacking_top5_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("Top5重要基因已保存至：", stacking_top5_path, "\n")

