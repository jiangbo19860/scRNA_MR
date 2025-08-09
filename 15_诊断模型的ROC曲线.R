# 输入与诺模图的输入一样：2个：ML模型输出的评分最高的基因文件3_outputs/20250807_importanceGene.RF.txt，合并之后的矩阵文件：1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt

# 加载必要的R包
library(caret)       # 用于机器学习模型训练
library(DALEX)       # 用于模型解释
library(ggplot2)     # 绘图包
library(randomForest) # 随机森林模型
library(kernlab)     # SVM 模型
library(xgboost)     # XGBoost模型
library(pROC)        # ROC曲线和AUC计算

# 设置随机种子以保证结果可重复
set.seed(123)

# 设置输入文件路径
inputFile = "merge.normalize.txt"  # 输入的标准化数据文件
geneFile = "importanceGene.RF.txt"  # 重要基因文件

# 设置工作目录
setwd("/Users/lijiangbo/scRNA_MR/4_references/Stroke_MR+机器学习/33诊断模型的ROC曲线")  # 设置工作目录

# 读取标准化数据文件
data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
row.names(data) = gsub("-", "_", row.names(data))  # 将行名中的 "-" 替换为 "_"

# 读取重要基因文件，并提取相关基因数据
geneRT = read.table(geneFile, header = TRUE, sep = "\t", check.names = FALSE)
data = data[as.vector(geneRT[, 1]), ]  # 根据基因列表筛选数据

# 提取样本信息并构建数据框
data = t(data)  # 转置数据，使样本为行，基因为列
group = gsub("(.*)\\_(.*)", "\\2", row.names(data))  # 提取样本分组信息
data = as.data.frame(data)
data$Type = group  # 添加样本类别信息

# 数据集划分为训练集和测试集
inTrain <- createDataPartition(y = data$Type, p = 0.7, list = FALSE)  # 按70%划分训练集
train <- data[inTrain, ]  # 训练集
test <- data[-inTrain, ]  # 测试集

# 根据基因文件选择模型进行训练
control = trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE)  # 设置5折交叉验证
if(geneFile == "importanceGene.RF.txt") {
  # 随机森林模型训练
  model = train(Type ~ ., data = train, method = 'rf', trControl = control)
} else if(geneFile == "importanceGene.SVM.txt") {
  # SVM模型训练
  model = train(Type ~ ., data = train, method = "svmRadial", prob.model = TRUE, trControl = control)
} else if(geneFile == "importanceGene.XGB.txt") {
  # XGBoost模型训练
  model = train(Type ~ ., data = train, method = "xgbDART", trControl = control)
} else if(geneFile == "importanceGene.GLM.txt") {
  # GLM模型训练
  model = train(Type ~ ., data = train, method = "glm", family = "binomial", trControl = control)
}

# 生成ROC曲线
yTest = ifelse(test$Type == "Control", 0, 1)  # 将类别转换为0和1
pred1 = predict(model, newdata = test, type = "prob")  # 预测测试集
roc1 = roc(yTest, as.numeric(pred1[, 2]))  # 计算ROC曲线
ci1 = ci.auc(roc1, method = "bootstrap")  # 计算AUC的95%置信区间
ciVec = as.numeric(ci1)  # 将置信区间转换为数值向量

# 保存ROC曲线到PDF文件
pdf(file = "ROC.pdf", width = 5, height = 5)
plot(roc1, print.auc = TRUE, legacy.axes = TRUE, main = "", col = "red")  # 绘制ROC曲线
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")  # 添加95% CI注释
dev.off()  # 关闭PDF输出
