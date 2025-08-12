rm(list = ls())
# 加载所需R包（新增sva包用于批次校正，tidyr用于数据重塑）
pacman::p_load(
  caret, randomForest, kernlab, xgboost, pROC, glmnet,
  DALEX, DALEXtra, ggplot2, here, dplyr, tibble, gridExtra,
  sva, tidyr  # 用于ComBat批次校正和数据处理
)

# 基础设置
set.seed(123)
today_date <- format(Sys.Date(), "%Y%m%d")
cat("分析日期：", today_date, "\n\n")

# 输入路径
inputFile <- here("1_data/GEO/GSE183795/GSE183795_normalized_matrix_final.txt")
diff_gene_file <- here("3_outputs/20250810_original_diffGeneExp_p05_26genes.txt")

# 输出路径（用列表统一管理）
out_paths <- list(
  sig_genes = here("3_outputs", paste0(today_date, "_extracted_sig_genes.txt")),
  residual = here("3_outputs", paste0(today_date, "_ML_residual.pdf")),  # 残差累积分布图
  boxplot = here("3_outputs", paste0(today_date, "_ML_boxplot.pdf")),    # 残差箱线图
  roc = here("3_outputs", paste0(today_date, "_ML_ROC.pdf")),
  importance = here("3_outputs", paste0(today_date, "_ML_importance.pdf")),  # 特征重要性图
  rf_genes = here("3_outputs", paste0(today_date, "_importanceGene.RF.txt")),
  svm_genes = here("3_outputs", paste0(today_date, "_importanceGene.SVM.txt")),
  xgb_genes = here("3_outputs", paste0(today_date, "_importanceGene.XGB.txt")),
  glm_genes = here("3_outputs", paste0(today_date, "_importanceGene.GLM.txt")),
  stable_genes = here("3_outputs", paste0(today_date, "_stable_importance_genes.txt")),
  stacking_pred = here("3_outputs", paste0(today_date, "_stacking_predictions.txt")),
  stacking_roc = here("3_outputs", paste0(today_date, "_stacking_ROC.pdf")),
  stacking_genes = here("3_outputs", paste0(today_date, "_importanceGene.Stacking.txt"))
)


### 1. 数据输入与预处理 ###
# 读取差异基因
diff_genes <- read.delim(diff_gene_file, header = TRUE, sep = "\t", check.names = FALSE)
sig_genes <- diff_genes$ID
cat("显著基因总数：", length(sig_genes), "（前6个：", paste(head(sig_genes), collapse = ", "), "）\n")
write.table(data.frame(Gene = sig_genes), out_paths$sig_genes,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 读取表达数据并去重
expr_data <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
colnames(expr_data)[1] <- "gene_name"
duplicate_genes <- unique(expr_data$gene_name[duplicated(expr_data$gene_name)])
cat("重复基因数量：", length(duplicate_genes), "\n")

# 合并重复基因（取均值）
expr_unique <- if (length(duplicate_genes) > 0) {
  bind_rows(
    expr_data %>% filter(gene_name %in% duplicate_genes) %>% group_by(gene_name) %>% summarise(across(everything(), \(x) mean(x, na.rm = TRUE))),
    expr_data %>% filter(!gene_name %in% duplicate_genes)
  )
} else {
  expr_data
}
rownames(expr_unique) <- expr_unique$gene_name
expr_unique <- expr_unique[, -1, drop = FALSE]
cat("去重后表达数据维度：", nrow(expr_unique), "×", ncol(expr_unique), "\n")

# 标准化基因名
rownames(expr_unique) <- gsub("-|\\.", "_", rownames(expr_unique)) %>% toupper()
sig_genes <- gsub("-|\\.", "_", sig_genes) %>% toupper()

# 筛选样本和重叠基因
expr_filtered <- expr_unique[, !grepl("NA", colnames(expr_unique)), drop = FALSE]
cat("筛选后样本数：", ncol(expr_filtered), "\n")
existing_genes <- intersect(sig_genes, rownames(expr_filtered))
if (length(existing_genes) == 0) stop("无匹配的显著基因！")
expr_sig <- expr_filtered[existing_genes, , drop = FALSE]
cat("用于建模的显著基因数量：", nrow(expr_sig), "\n")


### 2. 批次效应校正 ###
# 定义批次划分函数（根据样本名称前缀）
extract_batch <- function(sample_name) {
  if (grepl("^X\\d+", sample_name)) {
    return(gsub("^(X\\d+).*", "\\1", sample_name))
  } else if (grepl("^E\\d+", sample_name)) {
    return(gsub("^(E\\d+).*", "\\1", sample_name))
  } else if (grepl("^Hussain", sample_name)) {
    return("Hussain")
  } else {
    return("Other")
  }
}

# 提取批次信息和样本分组
sample_names <- colnames(expr_sig)
batch <- sapply(sample_names, extract_batch)
batch <- factor(batch)
sample_groups <- factor(gsub(".+_([^_]+)$", "\\1", sample_names), levels = c("Control", "Tumor"))

cat("批次数量：", length(levels(batch)), "\n")
cat("主要批次分布：\n")
print(head(sort(table(batch), decreasing = TRUE), 5))

# 执行ComBat批次校正
cat("执行批次效应校正...\n")
mod <- model.matrix(~sample_groups)
combat_expr <- ComBat(dat = expr_sig, batch = batch, mod = mod)
expr_corrected <- combat_expr


### 3. 异常样本检测与移除 ###
cat("Detecting and removing outliers...\n")
pca <- prcomp(t(expr_corrected), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])
distances <- dist(pca_df) %>% as.matrix() %>% rowMeans()
outlier_threshold <- quantile(distances, 0.95)
normal_samples <- names(distances)[distances < outlier_threshold]

# 更新表达矩阵、样本名和样本分组
expr_filtered <- expr_corrected[, normal_samples, drop = FALSE]
sample_names <- sample_names[sample_names %in% normal_samples]
cat("异常样本移除后剩余样本数：", length(sample_names), "\n")

# 重新提取样本分组
sample_groups <- sapply(sample_names, function(name) {
  if (grepl("_Control$", name)) return("Control")
  if (grepl("_Tumor$", name)) return("Tumor")
  warning(paste("样本名", name, "格式不符合预期，无法提取组别"))
  return(NA)
})
sample_groups <- factor(sample_groups, levels = c("Control", "Tumor"))

# 划分训练集和测试集
inTrain <- createDataPartition(y = sample_groups, p = 0.7, list = FALSE)
train_samples <- sample_names[inTrain]
test_samples <- sample_names[-inTrain]
cat("Training set samples: ", length(train_samples),
    "; Test set samples: ", length(test_samples), "\n", sep = "")


### 4. 数据转换与特征预处理 ###
# 转置矩阵并提取组别
data_samples <- as.data.frame(t(expr_filtered))
data_samples$Type <- sample_groups
cat("样本组别分布：", paste(names(table(data_samples$Type)), table(data_samples$Type), sep = ":", collapse = " | "), "\n")

# 特征预处理（去低方差和高相关特征）
feature_cols <- setdiff(colnames(data_samples), "Type")
feature_vars <- apply(data_samples[, feature_cols], 2, var, na.rm = TRUE)
zero_var <- sum(feature_vars == 0)
if (zero_var > 0) {
  data_samples <- data_samples[, c(names(feature_vars)[feature_vars > 0], "Type"), drop = FALSE]
  feature_cols <- setdiff(colnames(data_samples), "Type")
  cat("删除", zero_var, "个方差为0的特征，剩余：", length(feature_cols), "\n")
}
if (length(feature_cols) == 0) stop("所有特征方差为0，无法建模！")

# 去除高相关特征（相关系数>0.7）
if (length(feature_cols) > 1) {
  cor_mat <- cor(data_samples[, feature_cols])
  high_cor <- findCorrelation(cor_mat, cutoff = 0.7)
  if (length(high_cor) > 0) {
    data_samples <- data_samples[, c(feature_cols[-high_cor], "Type"), drop = FALSE]
    feature_cols <- feature_cols[-high_cor]
    cat("删除", length(high_cor), "个高相关特征，剩余：", length(feature_cols), "\n")
  }
}


### 5. 数据集划分与模型训练 ###
# 使用全局划分的样本索引
train_data <- data_samples[train_samples, , drop = FALSE]
test_data <- data_samples[test_samples, , drop = FALSE]
cat("训练集：", nrow(train_data), "；测试集：", nrow(test_data), "\n")

# 交叉验证配置
cv_control <- trainControl(
  method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE,
  summaryFunction = twoClassSummary, sampling = "up", verboseIter = FALSE
)

# 训练4个基模型（用列表存储模型）
models <- list()

# 随机森林
cat("训练随机森林模型...\n")
models[["RF"]] <- train(Type ~ ., data = train_data, method = "rf", trControl = cv_control,
                        metric = "ROC", ntree = 500)

# SVM
cat("训练SVM模型...\n")
models[["SVM"]] <- train(Type ~ ., data = train_data, method = "svmRadial", trControl = cv_control,
                         metric = "ROC", tuneLength = 5)

# XGBoost
cat("训练XGBoost模型...\n")
models[["XGB"]] <- train(Type ~ ., data = train_data, method = "xgbTree", trControl = cv_control,
                         metric = "ROC", tuneGrid = expand.grid(
                           nrounds = 100, eta = 0.3, max_depth = 3, gamma = 0,
                           colsample_bytree = 0.7, min_child_weight = 1, subsample = 0.8
                         ))

# GLM（单独处理特征）
cat("训练GLM模型...\n")
cor_mat <- cor(train_data[, feature_cols])
high_cor_glm <- findCorrelation(cor_mat, cutoff = 0.8)
glm_features <- if (length(high_cor_glm) > 0) feature_cols[-high_cor_glm] else feature_cols
max_glm_feat <- min(floor(nrow(train_data)/2), length(glm_features))
if (max_glm_feat < length(glm_features)) {
  anova_p <- summary(aov(as.formula(paste("Type ~", paste(glm_features, collapse = " + "))), data = train_data))[[1]][, "Pr(>F)"]
  glm_features <- names(sort(anova_p)[1:max_glm_feat])
}
cat("GLM使用特征数：", length(glm_features), "\n")
models[["GLM"]] <- train(Type ~ ., data = train_data[, c(glm_features, "Type")], method = "glmnet",
                         family = "binomial", trControl = cv_control, metric = "ROC",
                         tuneGrid = expand.grid(alpha = c(0, 0.5, 1), lambda = seq(0.001, 0.1, length.out = 10)))


### 6. 模型评估 ###
y_test <- as.numeric(test_data$Type == "Tumor")
pred_list <- lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  newdata <- if (model_name == "GLM") test_data[, c(glm_features, "Type")] else test_data
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
})
names(pred_list) <- names(models)

# 计算ROC和AUC
roc_list <- lapply(names(pred_list), function(model_name) {
  roc(y_test, pred_list[[model_name]])
})
names(roc_list) <- names(models)

# 输出评估结果
cat("\n模型AUC：\n")
sapply(names(roc_list), function(model_name) {
  cat(model_name, "：", roc_list[[model_name]]$auc, "\n")
})

# AUC置信区间和Youden指数
cat("\nAUC 95%置信区间：\n")
youden_list <- list()
sapply(names(roc_list), function(model_name) {
  roc_obj <- roc_list[[model_name]]
  ci <- ci.auc(roc_obj)
  cat(model_name, "：", roc_obj$auc, " [", round(ci[1], 3), ",", round(ci[3], 3), "]\n")
  youden_list[[model_name]] <<- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
})

cat("\n最优阈值性能：\n")
sapply(names(roc_list), function(model_name) {
  roc_obj <- roc_list[[model_name]]
  youden <- youden_list[[model_name]]
  cat(model_name, "：灵敏度=", round(roc_obj$sensitivities[youden], 3),
      "，特异度=", round(roc_obj$specificities[youden], 3), "\n")
})


### 7. 残差分析（绘制ML_boxplot图） ###
cat("进行残差分析...\n")

# 创建模型解释器（用于残差计算）
data_explain <- train_data[, feature_cols, drop = FALSE]
y_explain <- as.numeric(train_data$Type == "Tumor")
explainers <- lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  predict_fun <- function(m, newdata) {
    predict(m, newdata = newdata, type = "prob")[, "Tumor"]
  }
  data <- if (model_name == "GLM") {
    data_explain[, intersect(colnames(data_explain), glm_features), drop = FALSE]
  } else {
    data_explain
  }
  DALEX::explain(
    model = model,
    data = data,
    y = y_explain,
    predict_function = predict_fun,
    label = model_name,
    type = "classification"
  )
})
names(explainers) <- names(models)

# 计算模型性能（残差等指标）
mp_list <- lapply(explainers, model_performance)

# 输出1：残差累积分布图（保存到residual路径）
pdf(file = out_paths$residual, width = 6, height = 6)
p_residual_cum <- plot(mp_list[["RF"]], mp_list[["SVM"]], mp_list[["XGB"]], mp_list[["GLM"]])
print(p_residual_cum)
dev.off()
cat("残差累积分布图已保存至：", out_paths$residual, "\n")

# 输出2：残差箱线图（保存到boxplot路径，即ML_boxplot）
pdf(file = out_paths$boxplot, width = 6, height = 6)
p_residual_box <- plot(mp_list[["RF"]], mp_list[["SVM"]], mp_list[["XGB"]], mp_list[["GLM"]], geom = "boxplot")
print(p_residual_box)
dev.off()
cat("残差箱线图（ML_boxplot）已保存至：", out_paths$boxplot, "\n")


### 8. 绘制基模型的ROC曲线 ###
cat("绘制基模型ROC曲线并保存...\n")
pdf(out_paths$roc, width = 6, height = 6)
plot(roc_list[["RF"]], col = "red", lwd = 2, legacy.axes = TRUE,
     main = "ROC Curves of Base Models",
     xlab = "1 - Specificity", ylab = "Sensitivity")
plot(roc_list[["SVM"]], col = "blue", lwd = 2, add = TRUE)
plot(roc_list[["XGB"]], col = "green", lwd = 2, add = TRUE)
plot(roc_list[["GLM"]], col = "purple", lwd = 2, add = TRUE)
legend("bottomright",
       legend = paste0(names(roc_list), " (AUC: ", round(sapply(roc_list, function(x) x$auc), 3), ")"),
       col = c("red", "blue", "green", "purple"),
       lwd = 2, bty = "n", cex = 0.8)
dev.off()
cat("基模型ROC曲线已保存至：", out_paths$roc, "\n")


### 9. 特征重要性分析（绘制ML_importance图） ###
cat("进行特征重要性分析...\n")

# 提取各模型的特征重要性（Top25）
extract_importance <- function(explainer, top_n = 25) {
  imp <- variable_importance(explainer, loss_function = loss_root_mean_square)
  imp <- imp[imp$variable != "_baseline_", ]
  if (nrow(imp) > 0) imp[order(-abs(imp$dropout_loss)), ][1:top_n, ] else data.frame()
}

summarize_importance <- function(raw_imp) {
  if (nrow(raw_imp) == 0) return(data.frame())
  raw_imp %>% group_by(variable, label) %>%
    summarise(mean_dropout_loss = mean(dropout_loss), n_iterations = n(), .groups = "drop") %>%
    arrange(-mean_dropout_loss)
}

# 提取并汇总各模型重要性
imp_raw_list <- lapply(explainers, extract_importance, top_n = 25)
imp_summary_list <- lapply(imp_raw_list, summarize_importance)

# 保存各模型重要基因
imp_paths <- list(RF = out_paths$rf_genes, SVM = out_paths$svm_genes,
                  XGB = out_paths$xgb_genes, GLM = out_paths$glm_genes)
lapply(names(imp_summary_list), function(model_name) {
  imp_summary <- imp_summary_list[[model_name]]
  if (nrow(imp_summary) > 0) {
    write.table(imp_summary[, c("variable", "mean_dropout_loss", "n_iterations", "label")],
                imp_paths[[model_name]], sep = "\t", quote = FALSE, row.names = FALSE)
    cat(model_name, "重要基因已保存至：", imp_paths[[model_name]], "\n")
  }
})

# 绘制特征重要性图（ML_importance）
cat("绘制特征重要性图...\n")
pdf(file = out_paths$importance, width = 12, height = 10)

# 为每个模型绘制Top10重要特征条形图
plots <- lapply(names(imp_summary_list), function(model_name) {
  imp_df <- imp_summary_list[[model_name]]
  if (nrow(imp_df) == 0) return(NULL)

  # 取Top10特征
  top10 <- imp_df %>% head(10) %>%
    mutate(variable = factor(variable, levels = rev(variable)))  # 倒序排列

  # 绘制条形图
  ggplot(top10, aes(x = variable, y = mean_dropout_loss)) +
    geom_bar(stat = "identity", fill = ifelse(model_name == "RF", "#FF6B6B",
                                              ifelse(model_name == "SVM", "#4ECDC4",
                                                     ifelse(model_name == "XGB", "#45B7D1", "#FFA07A")))) +
    coord_flip() +
    labs(title = paste(model_name),
         x = "", y = "Mean dropout loss") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
})

# 组合图形（过滤空图）
plots <- plots[!sapply(plots, is.null)]
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 2)  # 2列布局
}
dev.off()
cat("特征重要性图（ML_importance）已保存至：", out_paths$importance, "\n")


### 10. 筛选稳定基因 ###
read_importance <- function(file_path, model_name) {
  if (file.exists(file_path)) {
    df <- read.delim(file_path, sep = "\t", stringsAsFactors = FALSE)
    if ("variable" %in% colnames(df) && "mean_dropout_loss" %in% colnames(df)) {
      return(data.frame(variable = df$variable, mean_dropout_loss = df$mean_dropout_loss, model = model_name))
    }
  }
  data.frame()
}

all_imp <- bind_rows(lapply(names(imp_paths), function(model_name) {
  read_importance(imp_paths[[model_name]], model_name)
}))

if (nrow(all_imp) > 0) {
  gene_freq <- all_imp %>% group_by(variable) %>%
    summarise(occurrence = n_distinct(model), mean_importance = mean(mean_dropout_loss), .groups = "drop") %>%
    arrange(desc(occurrence), desc(mean_importance))

  # 确保至少5个稳定基因
  target_count <- 5
  stable_genes <- gene_freq %>% filter(occurrence == 4)
  if (nrow(stable_genes) < target_count) {
    stable_genes <- bind_rows(stable_genes, gene_freq %>% filter(occurrence == 3) %>% head(target_count - nrow(stable_genes)))
  }
  if (nrow(stable_genes) < target_count) {
    stable_genes <- bind_rows(stable_genes, gene_freq %>% filter(occurrence == 2) %>% head(target_count - nrow(stable_genes)))
  }
  if (nrow(stable_genes) < target_count) {
    stable_genes <- bind_rows(stable_genes, gene_freq %>% filter(occurrence == 1) %>% head(target_count - nrow(stable_genes)))
  }

  write.table(stable_genes, out_paths$stable_genes, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("稳定重要基因已保存至：", out_paths$stable_genes, "\n")
}


### 11. Stacking集成模型 ###
cat("\n训练Stacking集成模型...\n")
# 提取基模型预测概率
train_pred <- lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  newdata <- if (model_name == "GLM") train_data[, c(glm_features, "Type")] else train_data
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
})
names(train_pred) <- names(models)
train_stacking <- data.frame(train_pred, Type = train_data$Type)

test_pred <- lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  newdata <- if (model_name == "GLM") test_data[, c(glm_features, "Type")] else test_data
  predict(model, newdata = newdata, type = "prob")[, "Tumor"]
})
names(test_pred) <- names(models)
test_stacking <- data.frame(test_pred)

# 训练元模型
meta_model <- train(Type ~ ., data = train_stacking, method = "glmnet", family = "binomial",
                    trControl = cv_control, metric = "ROC",
                    tuneGrid = expand.grid(alpha = c(0, 0.5, 1), lambda = seq(0.001, 0.1, length.out = 10)))

# 评估Stacking模型
meta_pred_prob <- predict(meta_model, newdata = test_stacking, type = "prob")[, "Tumor"]
roc_stacking <- roc(y_test, meta_pred_prob)
cat("Stacking模型AUC：", roc_stacking$auc, "\n")

# 保存Stacking结果
write.table(data.frame(
  Sample = rownames(test_data), Actual = test_data$Type,
  Stacking_Prob = meta_pred_prob,
  Stacking_Pred = ifelse(meta_pred_prob >= roc_stacking$thresholds[which.max(roc_stacking$sensitivities + roc_stacking$specificities - 1)], "Tumor", "Control")
), out_paths$stacking_pred, sep = "\t", quote = FALSE, row.names = FALSE)

# 绘制含Stacking的ROC曲线
pdf(out_paths$stacking_roc, width = 6, height = 6)
plot(roc_list[["RF"]], col = "red", lwd = 2, legacy.axes = TRUE, main = "ROC Curves (Including Stacking)")
plot(roc_list[["SVM"]], col = "blue", lwd = 2, add = TRUE)
plot(roc_list[["XGB"]], col = "green", lwd = 2, add = TRUE)
plot(roc_list[["GLM"]], col = "purple", lwd = 2, add = TRUE)
plot(roc_stacking, col = "black", lwd = 2, lty = 2, add = TRUE)
legend("bottomright",
       legend = c(paste0(names(roc_list), " (AUC: ", round(sapply(roc_list, function(x) x$auc), 3), ")"),
                  paste0("Stacking (AUC: ", round(roc_stacking$auc, 3), ")")),
       col = c("red", "blue", "green", "purple", "black"),
       lwd = 2, lty = c(1, 1, 1, 1, 2), bty = "n", cex = 0.8)
dev.off()

# Stacking模型重要基因
base_importance <- bind_rows(lapply(names(imp_paths), function(model_name) {
  df <- read.delim(imp_paths[[model_name]], sep = "\t", stringsAsFactors = FALSE)
  if (nrow(df) > 0) df %>% select(variable, mean_dropout_loss) %>% mutate(model = model_name) else data.frame()
}))

meta_weights <- data.frame(
  model = names(models),
  weight = abs(coef(meta_model$finalModel, meta_model$bestTune$lambda)[-1, 1])
)
meta_weights$weight <- meta_weights$weight / sum(meta_weights$weight)

stacking_imp <- base_importance %>% left_join(meta_weights, by = "model") %>%
  mutate(weighted_importance = mean_dropout_loss * weight) %>%
  group_by(variable) %>% summarise(total_importance = mean(weighted_importance, na.rm = TRUE)) %>%
  arrange(desc(total_importance)) %>% head(5)

write.table(stacking_imp, out_paths$stacking_genes, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Stacking模型Top5重要基因已保存至：", out_paths$stacking_genes, "\n")

cat("\n===== 所有输出文件保存路径 =====", "\n")
for (file_label in names(out_paths)) {
  if (file.exists(out_paths[[file_label]])) {
    cat(sprintf("%-15s: %s\n", file_label, out_paths[[file_label]]))
  } else {
    cat(sprintf("%-15s: %s（文件未生成）\n", file_label, out_paths[[file_label]]))
  }
}
