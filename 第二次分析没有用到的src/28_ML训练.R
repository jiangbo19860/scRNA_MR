# 目标变量是Type(Control or Tumor)? 4种模型：随机森林、SVM、XGBoost、GLM均是用于二分类任务
# 基于基因表达数据和显著基因列表，构建并评估多种机器学习模型，以实现样本组别的分类（如疾病 vs 对照），同时分析关键基因的重要性。代码主要分为数据预处理、模型训练、模型评估和结果输出四个环节。
# 2个输入文件：合并的标准化基因表达数据1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt和显著意义的基因列表（3_outputs/汇总有显著意义的MR结果/gene.txt - geneRT）。
# 三、输出文件
# 评估图表
# ML_residual.pdf：模型残差的累积分布图，反映不同模型的预测误差分布。
# ML_boxplot.pdf：模型残差的箱线图，直观比较不同模型的误差离散程度。
# ML_ROC.pdf：ROC 曲线及 AUC 值，比较 4 种模型的分类性能（曲线越靠近左上角，AUC 越大，性能越好）。
# ML_importance.pdf：特征重要性图，展示各模型中对分类贡献最大的基因及其重要性分值。
# 重要基因列表
# importanceGene.RF.txt：随机森林模型中前 5 个最重要的基因。
# importanceGene.SVM.txt：支持向量机模型中前 5 个最重要的基因。
# importanceGene.XGB.txt：XGBoost 模型中前 5 个最重要的基因。
# importanceGene.GLM.txt：广义线性模型中前 5 个最重要的基因。

rm(list = ls())  # 清理工作环境
pacman::p_load(
  caret,      # 用于机器学习的训练控制和数据分割
  DALEX,      # 用于模型解释的包
  ggplot2,    # 用于绘图
  randomForest, # 随机森林算法
  kernlab,    # 支持向量机算法
  xgboost,    # XGBoost 算法
  pROC,        # 用于绘制和评估ROC曲线
  here
)

# 设置随机种子以保证结果可重复
set.seed(123)

# 设置输入文件路径
inputFile = here("1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt")  # 标准化的合并数据
geneFile = here("3_outputs/20250804_sig_genes_44.txt")  # 基因列表文件

# 读取基因表达数据文件（行：基因，列：样本）
data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 筛选出列名中不包含"NA"的列（保留有效样本列）
data <- data[, !grepl("NA", colnames(data)), drop = FALSE]
cat("清理后原始数据的样本列数：", ncol(data), "\n")

# 读取基因列表文件，并提取所需的基因表达数据
geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data_gene = data[as.vector(geneRT[, 1]), , drop = FALSE]  # 提取目标基因（行：基因，列：样本）

# 删除因基因名不匹配导致的全NA基因行（无效基因）
data_gene <- data_gene[!apply(is.na(data_gene), 1, all), , drop = FALSE]
cat("删除全NA基因后，剩余有效基因数量：", nrow(data_gene), "\n")
if (nrow(data_gene) == 0) stop("无有效基因，请检查gene.txt中的基因名是否与原始数据匹配！")

# 处理基因名中的特殊字符
row.names(data_gene) = gsub("-", "_", row.names(data_gene))

# 转置数据：行→样本，列→基因（便于后续按样本处理）
data_transposed = t(data_gene)  # 转置后：行是样本，列是基因
data_transposed = as.data.frame(data_transposed)

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

# 为清理后的的样本添加组别信息（从样本ID提取Control/Tumor）
data_clean$Type <- gsub("(.*)\\_(.*)", "\\2", rownames(data_clean))  # 提取组别标签
cat("清理后样本的组别分布：", table(data_clean$Type), "\n")  # 确认组别信息正确

# 检查并删除方差为0的特征（无区分度的基因）
feature_vars <- apply(data_clean[, -ncol(data_clean)], 2, var)
zero_var_features <- sum(feature_vars == 0)
cat("方差为0的特征数量：", zero_var_features, "\n")

if (zero_var_features > 0) {
  valid_features <- names(feature_vars)[feature_vars > 0]
  data_clean <- data_clean[, c(valid_features, "Type")]
  cat("删除无区分度特征后，剩余有效特征数量：", length(valid_features), "\n")
  if (length(valid_features) == 0) stop("所有特征均无区分度，无法训练模型！")
}

colnames(data_clean)  # 查看清理后的数据列名（应包含基因名和Type列）

# 划分独立的 “训练集” 和 “测试集”，用于最终评估模型的泛化能力。
inTrain <- createDataPartition(y = data_clean$Type, p = 0.7, list = FALSE)  # 使用清理后的数据
train <- data_clean[inTrain, ]  # 训练集
test <- data_clean[-inTrain, ]  # 测试集
cat("训练集样本数：", nrow(train), "，测试集样本数：", nrow(test), "\n")
cat("训练集组别分布：", table(train$Type), "\n")

# 训练集内部用留一法（LOOCV）进行交叉验证，用于模型调优（如选择最优参数）。留一法的本质是 “样本量 = 折数”。训练集24个样本量等价于24折交叉验证（样本量=24时），每次用n-1个样本训练，1个样本验证
control = trainControl(
  method = "LOOCV",  # 改为留一法交叉验证
  savePredictions = TRUE,
  classProbs = TRUE  # 保留类别概率，用于后续ROC分析
)

# 随机森林模型训练
mod_rf = train(Type ~ ., data = train, method = 'rf', trControl = control)  # mod_rf是模型对象的名称，用于存储训练完成后的随机森林模型。Type：目标变量（需要预测的变量），即样本的组别（Control 或 Tumor）。~：表示 “由…… 预测”。.：表示 “除目标变量外的所有其他变量”，即 train 数据框中除 Type 列之外的所有列（这里是基因表达特征）。整体含义：用所有基因特征预测样本的组别（Type）。'rf' 表示使用随机森林（Random Forest）算法。trControl = control传递训练控制参数，control 是之前定义的交叉验证策略（留一法）。

# 支持向量机模型训练
mod_svm = train(Type ~ ., data = train, method = "svmRadial", prob.model = TRUE, trControl = control)  # 训练SVM模型


# XGBoost模型训练_如果样本量太小，要减少特征数，最好是样本数：特征数=5~10：1 ----------------------------
# 计算特征与类别的相关性（以ANOVA为例）
feature_anova <- apply(train[, -ncol(train), drop = FALSE], 2, function(x) { # drop=FALSE是保持原始数据的维度结构，避免因子集提取导致维度意外简化。
  # 用tryCatch捕获异常，避免单个特征错误导致整体为NULL
  tryCatch({
    model <- aov(x ~ train$Type)
    summary(model)[[1]]$`Pr(>F)`[1]  # 更稳定的P值提取方式
  }, error = function(e) {
    NA  # 异常特征标记为NA，而非导致整体为NULL
  })
})
feature_anova
# class(train_sub$Type)  # 应为"factor"，而非"character"
# # 将Type转换为因子
# train_sub$Type <- as.factor(train_sub$Type)
# test_sub$Type <- as.factor(test_sub$Type)
# class(train_sub$Type)  # 确认输出"factor"

# 保留P值最小的前5个特征（最具区分度）
top_features <- names(sort(feature_anova)[1:min(5, length(feature_anova))])
train_sub <- train[, c(top_features, "Type")]  # 子集训练集
test_sub <- test[, c(top_features, "Type")]    # 子集测试集

# 检查train_sub的类型
class(train_sub)  # 应输出"data.frame"

# 若不是数据框，强制转换
if (!inherits(train_sub, "data.frame")) {
  train_sub <- as.data.frame(train_sub)
  test_sub <- as.data.frame(test_sub)
}
# 查看feature_anova是否有效（应包含特征名和对应的P值）
head(feature_anova)
# 查看top_features是否为空（正常应包含多个特征名）
top_features
# 确认特征筛选后的数据是否完整
head(train_sub)  # 查看前几行
ncol(train_sub)  # 应比原train多1列（包含Type）

# 查看train_sub的类别分布
table(train_sub$Type)

# 若某类样本数≤1，留一法必然出错（需调整数据划分）
if (min(table(train_sub$Type)) <= 1) {
  # 降低训练集比例（如60%），确保两类样本数≥2
  inTrain <- createDataPartition(y = data_clean$Type, p = 0.6, list = FALSE)
  train <- data_clean[inTrain, ]
  train_sub <- train[, c(top_features, "Type")]  # 重新生成train_sub
  cat("调整后训练集类别分布：", table(train_sub$Type), "\n")
}
# 确认Type是两水平因子（无多余水平）
levels(train_sub$Type)  # 应输出"Control"和"Tumor"

# 检查是否有异常值（极端值可能导致模型无法收敛）
summary(train_sub[, feature_vars])  # 若某特征存在极大/极小值，可适当截断或标准化
# 1. 将Type转换为因子型，并指定水平为"Control"和"Tumor"
train_sub$Type <- factor(train_sub$Type, levels = c("Control", "Tumor"))
test_sub$Type <- factor(test_sub$Type, levels = c("Control", "Tumor"))  # 测试集同步转换

# 2. 验证转换结果
levels(train_sub$Type)  # 应输出"Control" "Tumor"
class(train_sub$Type)   # 应输出"factor"
# 1. 重新定义特征列名（排除Type列）
feature_vars <- setdiff(colnames(train_sub), "Type")  # 自动提取所有特征列名

# 2. 验证特征列是否存在
feature_vars  # 应输出5个特征名："NFATC2IP" "CLEC16A" "CSDC2" "HMCN1" "PLCE1"
all(feature_vars %in% colnames(train_sub))  # 应返回TRUE，说明列名匹配

# 3. 再次查看特征的统计摘要
summary(train_sub[, feature_vars])  # 此时应正常输出各特征的统计量（最小值、最大值等）

# 重新定义交叉验证策略（推荐留一法，适配小样本）
# control <- trainControl(
#   method = "LOOCV",  # 留一法交叉验证
#   classProbs = TRUE,  # 保留类别概率
#   summaryFunction = twoClassSummary  # 二分类专用评估函数（避免Accuracy缺失）
# )

# 定义参数网格（包含所有XGBoost参数，包括nrounds）
xgb_grid <- expand.grid(
  nrounds = 150,          # 表示最终模型由 150 棵决策树组成。
  eta = 0.4,              # 学习率
  max_depth = 3,          # 树深度，是指每棵决策树从根节点到叶节点的最长路径包含的节点数（）
  gamma = 0,              # 分裂门槛
  colsample_bytree = c(0.6, 0.8),  # 列采样
  min_child_weight = 1,   # 子节点最小权重
  subsample = c(0.5, 0.75, 1.0)    # 样本抽样参数，用于控制每棵树训练时随机抽取的样本比例，目的是减少过拟合并提高模型泛化能力。
)

# 训练模型时，只通过tuneGrid传递参数，不单独指定nrounds
mod_xgb <- train(
  Type ~ .,
  data = train_sub,
  method = "xgbTree",
  trControl = control,
  metric = "ROC",
  tuneGrid = xgb_grid  # 所有参数通过网格传递，包括nrounds
  # 注意：此处不要额外写nrounds = 150，否则会重复
)

# 广义线性模型（GLM）训练
mod_glm = train(Type ~ ., data = train, method = "glm", family = "binomial", trControl = control)  # 训练GLM模型

# 定义预测函数（用于模型解释）
p_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, 2]
}
yTest = ifelse(test$Type == "Control", 0, 1)  # 将测试集组别转换为0（Control）和1（Tumor）

# 随机森林模型的模型解释
explainer_rf = explain(mod_rf, label = "RF", data = test, y = yTest, predict_function = p_fun, verbose = FALSE)
mp_rf = model_performance(explainer_rf)

# 支持向量机模型的模型解释
explainer_svm = explain(mod_svm, label = "SVM", data = test, y = yTest, predict_function = p_fun, verbose = FALSE)
mp_svm = model_performance(explainer_svm)

# XGBoost模型的模型解释
explainer_xgb = explain(mod_xgb, label = "XGB", data = test, y = yTest, predict_function = p_fun, verbose = FALSE)
mp_xgb = model_performance(explainer_xgb)

# GLM模型的模型解释
explainer_glm = explain(mod_glm, label = "GLM", data = test, y = yTest, predict_function = p_fun, verbose = FALSE)
mp_glm = model_performance(explainer_glm)

# 生成模型残差的累积分布图
pdf(file = here("3_outputs/ML_residual.pdf"), width = 6, height = 6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

# 生成模型残差的箱线图
pdf(file = here("3_outputs/ML_boxplot.pdf"), width = 6, height = 6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()

# 绘制ROC曲线并比较模型性能
pred1 = predict(mod_rf, newdata = test, type = "prob")
pred2 = predict(mod_svm, newdata = test, type = "prob")
pred3 = predict(mod_xgb, newdata = test, type = "prob")
pred4 = predict(mod_glm, newdata = test, type = "prob")
# 为SVM的预测结果添加行名
rownames(pred2) <- rownames(test)

# 为XGBoost的预测结果添加行名
rownames(pred3) <- rownames(test)

roc1 = roc(yTest, as.numeric(pred1[, 2]))
roc2 = roc(yTest, as.numeric(pred2[, 2]))
roc3 = roc(yTest, as.numeric(pred3[, 2]))
roc4 = roc(yTest, as.numeric(pred4[, 2]))

pdf(file = here("3_outputs/ML_ROC.pdf"), width = 5, height = 5)
plot(roc1, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "red")
plot(roc2, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "blue", add = TRUE)
plot(roc3, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "green", add = TRUE)
plot(roc4, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "yellow", add = TRUE)
legend('bottomright',
       c(paste0('RF: ', sprintf("%.03f", roc1$auc)),
         paste0('SVM: ', sprintf("%.03f", roc2$auc)),
         paste0('XGB: ', sprintf("%.03f", roc3$auc)),
         paste0('GLM: ', sprintf("%.03f", roc4$auc))),
       col = c("red", "blue", "green", "yellow"), lwd = 2, bty = 'n')
dev.off()

# 提取重要特征并生成特征重要性图
importance_rf <- variable_importance(explainer_rf, loss_function = loss_root_mean_square)
importance_svm <- variable_importance(explainer_svm, loss_function = loss_root_mean_square)
importance_glm <- variable_importance(explainer_glm, loss_function = loss_root_mean_square)
importance_xgb <- variable_importance(explainer_xgb, loss_function = loss_root_mean_square)

# 绘制变量重要性图（调整索引确保取到有效特征）
n_features <- ncol(data_clean) - 1  # 总特征数（排除Type列）
plot_indices <- c(1, max(1, n_features - 8):n_features)  # 取前1个和后8个特征（避免索引越界）

pdf(file = here("3_outputs/ML_importance.pdf"), width = 7, height = 10)
plot(importance_rf[plot_indices, ],
     importance_svm[plot_indices, ],
     importance_xgb[plot_indices, ],
     importance_glm[plot_indices, ])
dev.off()

# 将最重要的基因写入文件（基于清理后的数据列数计算索引）
geneNum = 5  # 提取前5个重要基因
write.table(importance_rf[(n_features - geneNum + 1):n_features, ],
            file = here("3_outputs/importanceGene.RF.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_svm[(n_features - geneNum + 1):n_features, ],
            file = here("3_outputs/importanceGene.SVM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_xgb[(n_features - geneNum + 1):n_features, ],
            file = here("3_outputs/importanceGene.XGB.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_glm[(n_features - geneNum + 1):n_features, ],
            file = here("3_outputs/importanceGene.GLM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("所有分析完成！结果已保存至3_outputs文件夹。\n")
