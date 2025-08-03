# 安装必要的包（如果尚未安装）
#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("pROC")
#install.packages("xgboost")
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
geneFile = here("3_outputs/汇总有显著意义的MR结果/gene.txt")  # 基因列表文件

# 读取基因表达数据文件
data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# 读取原始文件的表头行（第一行）
header <- readLines(inputFile, n = 1)
# 按分隔符（\t）拆分表头，查看每个元素
strsplit(header, "\t")[[1]]
colnames(data)
# 检查原始数据（提取基因后）是否有缺失值
sum(is.na(data))  # 统计所有NA的数量
which(is.na(data), arr.ind = TRUE)  # 定位NA所在的行列（若数量较少）

# 筛选出列名中不包含"NA"的列（保留有效列）
data <- data[, !grepl("NA", colnames(data)), drop = FALSE]

# 验证结果（查看清理后的列名）
colnames(data)


# 读取基因列表文件，并提取所需的基因表达数据
geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data = data[as.vector(geneRT[, 1]), ]  # 提取基因表达数据
row.names(data) = gsub("-", "_", row.names(data))  # 将行名中的 "-" 替换为 "_"

# 提取样本的组别信息
data = t(data)  # 转置数据，将样本作为行
group = gsub("(.*)\\_(.*)", "\\2", row.names(data))  # 提取组别标签
data = as.data.frame(data)  # 将数据转换为数据框
data$Type = group  # 将组别信息添加到数据中

# 转置后、添加组别信息前，检查样本数据中的缺失值
data_transposed <- t(data)  # 转置后的数据（样本为行，基因为列）
cat("转置后样本数据中的缺失值数量：", sum(is.na(data_transposed)), "\n")

# 查看哪些样本（行）存在缺失值
na_samples <- which(apply(is.na(data_transposed), 1, any))
cat("存在缺失值的样本数量：", length(na_samples), "\n")

# 转置数据（样本为行，基因为列）
data_transposed <- t(data)
data_transposed <- as.data.frame(data_transposed)

# 查看含缺失值的样本ID
na_sample_ids <- rownames(data_transposed)[na_samples]
cat("含缺失值的样本ID：", paste(na_sample_ids, collapse = ", "), "\n")

# 删除含缺失值的样本（行）
data_clean <- data_transposed[!rownames(data_transposed) %in% na_sample_ids, , drop = FALSE]
cat("删除缺失值样本后，剩余样本数量：", nrow(data_clean), "\n")

# 验证缺失值是否已清除
cat("清理后样本数据中的缺失值数量：", sum(is.na(data_clean)), "\n")  # 应输出0

# 数据集划分为训练集和测试集
inTrain <- createDataPartition(y = data$Type, p = 0.7, list = FALSE)  # 按70%划分训练集
train <- data[inTrain, ]  # 训练集
test <- data[-inTrain, ]  # 测试集

# 随机森林模型训练
control = trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE)  # 设置5折交叉验证
mod_rf = train(Type ~ ., data = train, method = 'rf', trControl = control)  # 训练随机森林模型

# 支持向量机模型训练
mod_svm = train(Type ~ ., data = train, method = "svmRadial", prob.model = TRUE, trControl = control)  # 训练SVM模型

# XGBoost模型训练
mod_xgb = train(Type ~ ., data = train, method = "xgbDART", trControl = control)  # 训练XGBoost模型

# 广义线性模型（GLM）训练
mod_glm = train(Type ~ ., data = train, method = "glm", family = "binomial", trControl = control)  # 训练GLM模型

# 定义预测函数
p_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, 2]
}
yTest = ifelse(test$Type == "Control", 0, 1)  # 将测试集的类型转换为0和1

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

# 绘制ROC曲线
pred1 = predict(mod_rf, newdata = test, type = "prob")
pred2 = predict(mod_svm, newdata = test, type = "prob")
pred3 = predict(mod_xgb, newdata = test, type = "prob")
pred4 = predict(mod_glm, newdata = test, type = "prob")
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

# 绘制变量重要性图
pdf(file = here("3_outputs/ML_importance.pdf"), width = 7, height = 10)
plot(importance_rf[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_svm[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_xgb[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_glm[c(1, (ncol(data) - 8):(ncol(data) + 1)), ])
dev.off()

# 将最重要的基因写入文件
geneNum = 5  # 提取前5个重要基因
write.table(importance_rf[(ncol(data) - geneNum + 2):(ncol(data) + 1), ], file = here("3_outputs/importanceGene.RF.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_svm[(ncol(data) - geneNum + 2):(ncol(data) + 1), ], file = here("3_outputs/importanceGene.SVM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_xgb[(ncol(data) - geneNum + 2):(ncol(data) + 1), ], file = here("3_outputs/importanceGene.XGB.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_glm[(ncol(data) - geneNum + 2):(ncol(data) + 1), ], file = here("3_outputs/importanceGene.GLM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
