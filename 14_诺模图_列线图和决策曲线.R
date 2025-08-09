# 输入：2个：ML模型输出的评分最高的基因文件3_outputs/20250807_importanceGene.RF.txt，合并之后的矩阵文件：1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt。
# 需要修改的代码：lrmModel = lrm(Type ~FCER1G+	FHIT+	FARS2	+IL10RB	+CCDC88C, data = rt, x = TRUE, y = TRUE)  中的基因名要替换成自己的实际结果。dc = decision_curve(Type ~ FCER1G+	FHIT+	FARS2	+IL10RB	+CCDC88C,中的基因名也要修改。

# 安装必要的包（如果尚未安装）
#install.packages("rms")
#install.packages("rmda")

# 加载所需的R包
library(rms)  # 用于回归模型和Nomogram的包
library(rmda) # 用于决策曲线分析 (DCA) 的包

# 设置输入文件路径
inputFile = "merge.normalize.txt"  # 标准化数据文件
geneFile = "importanceGene.RF.txt"  # 重要基因文件


# 读取标准化数据文件
data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
row.names(data) = gsub("-", "_", row.names(data))  # 将行名中的 "-" 替换为 "_"

# 读取重要基因文件，并提取相关基因数据
geneRT = read.table(geneFile, header = TRUE, sep = "\t", check.names = FALSE)
data = data[as.vector(geneRT[, 1]), ]  # 根据基因列表筛选数据

# 提取样本信息并构建数据框
data = t(data)  # 转置数据，使样本为行，基因为列
group = gsub("(.*)\\_(.*)", "\\2", row.names(data))  # 提取样本分组信息
rt = cbind(as.data.frame(data), Type = group)  # 将基因数据和分组信息合并

# 数据准备，生成ddist对象
ddist = datadist(rt)  # 生成datadist对象，用于Nomogram建模
options(datadist = "ddist")  # 设置全局datadist选项

# 构建逻辑回归模型，并生成Nomogram
lrmModel = lrm(Type ~FCER1G+	FHIT+	FARS2	+IL10RB	+CCDC88C, data = rt, x = TRUE, y = TRUE)  # 构建逻辑回归模型
nomo = nomogram(lrmModel, fun = plogis, fun.at = c(0.0001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99), lp = FALSE, funlabel = "Risk of Disease")  # 生成Nomogram

# 保存Nomogram图像
pdf("Nomo.pdf", width = 8, height = 6)  # 打开PDF设备，保存图像
plot(nomo)  # 绘制Nomogram
dev.off()  # 关闭PDF设备

# 校准曲线的生成和绘制
cali = calibrate(lrmModel, method = "boot", B = 1000)  # 使用bootstrap方法进行模型校准
pdf("Calibration.pdf", width = 5.5, height = 5.5)  # 打开PDF设备，保存图像
plot(cali, xlab = "Predicted probability", ylab = "Actual probability", sub = FALSE)  # 绘制校准曲线
dev.off()  # 关闭PDF设备

# 决策曲线分析 (DCA)
rt$Type = ifelse(rt$Type == "Control", 0, 1)  # 将样本类别转为0和1
dc = decision_curve(Type ~ FCER1G+	FHIT+	FARS2	+IL10RB	+CCDC88C, data = rt,
                    family = binomial(link = 'logit'),
                    thresholds = seq(0, 1, by = 0.01),
                    confidence.intervals = 0.95)  # 生成决策曲线

# 绘制DCA曲线并保存
pdf(file = "DCA.pdf", width = 5.5, height = 5.5)  # 打开PDF设备，保存图像
plot_decision_curve(dc,
                    curve.names = "Model",
                    xlab = "Threshold probability",
                    cost.benefit.axis = TRUE,
                    col = "red",
                    confidence.intervals = FALSE,
                    standardize = FALSE)  # 绘制DCA曲线
dev.off()  # 关闭PDF设备
