# 2个输入：处理后的验证集GSE表达数据.txt, 和ML预测结果中最好的importanceGene.RF.txt。要新建一个gene.txt文档，把importanceGene.RF.txt中的gene名复制为一列。
# 输出是ROC曲线，
rm(list = ls())  # 清空工作空间
# 加载所需的包
library(glmnet)
library(pROC)

# 设置文件路径
expFile="GSE69528.normalize.txt"      # 表达数据文件
geneFile="gene.txt"         # 基因列表文件
setwd("/Users/lijiangbo/scRNA_MR/4_references/Stroke_MR+机器学习/34外部验证集诊断模型的ROC曲线")    # 设置工作目录

# 读取表达数据文件
rt = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# 提取样本的分组信息
y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y = ifelse(y == "Control", 0, 1)

# 读取基因列表文件
geneRT = read.table(geneFile, header=F, sep="\t", check.names=F)

# 循环遍历每个基因并生成ROC曲线
bioCol = rainbow(nrow(geneRT), s=0.9, v=0.9)    # 设置曲线的颜色
aucText = c()
k = 0
for (x in as.vector(geneRT[,1])) {
  k = k + 1
  # 生成ROC曲线
  roc1 = roc(y, as.numeric(rt[x,]))     # 获取ROC曲线的参数
  if (k == 1) {
    pdf(file="ROC.genes.pdf", width=5, height=4.75)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}
# 添加图例，显示每条ROC曲线对应的基因和AUC值
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

# 生成逻辑回归模型
rt = rt[as.vector(geneRT[,1]),]
rt = as.data.frame(t(rt))
logit = glm(y ~ ., family=binomial(link='logit'), data=rt)
pred = predict(logit, newx=rt)     # 获取模型的预测值

# 计算并绘制模型的ROC曲线
roc1 = roc(y, as.numeric(pred))      # 获取模型ROC曲线的参数
ci1 = ci.auc(roc1, method="bootstrap")     # 获取ROC曲线下面积的置信区间
ciVec = as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model")
# 填充ROC曲线下的面积
polygon(c(roc1$specificities, 1), c(roc1$sensitivities, 0), col=rgb(1, 0, 0, 0.2), border=NA)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.3f", ciVec[1]), "-", sprintf("%.3f", ciVec[3])), col="red")
dev.off()

# 随机扰动标签以降低AUC值
set.seed(123)  # 设置随机种子
y_noisy = sample(y)  # 随机打乱标签
logit_noisy = glm(y_noisy ~ ., family=binomial(link='logit'), data=rt)
pred_noisy = predict(logit_noisy, newx=rt)

# 计算并绘制扰动标签后的ROC曲线
roc1_noisy = roc(y_noisy, as.numeric(pred_noisy))
#pdf(file="ROC.model_noisy.pdf", width=5, height=4.75)
plot(roc1_noisy, print.auc=TRUE, col="blue", legacy.axes=T, main="改变临床标题")
# 填充ROC曲线下的面积
polygon(c(roc1_noisy$specificities, 1), c(roc1_noisy$sensitivities, 0), col=rgb(0, 0, 1, 0.2), border=NA)
dev.off()
