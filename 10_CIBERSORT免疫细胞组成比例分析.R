# CIBERSORT算法：解析复杂组织样本中的免疫细胞组成比例。CIBERSORT 是一种基于支持向量回归（Support Vector Regression, SVR）的计算方法，由美国斯坦福大学的研究团队于 2015 年提出，主要用于解析复杂组织样本中的免疫细胞组成比例。它通过利用已知的免疫细胞亚型特异性基因表达特征（称为 “特征矩阵”，signature matrix），反推混合样本（如肿瘤组织、血液等）中各种免疫细胞的相对丰度。
# input有2个：参考基因表达矩阵（Signature Matrix，通常为ref.txt）和混合样本的基因表达数据（Mixture File，通常为inputFile，行是基因（需与参考矩阵中的基因名匹配），列是样本（如患者样本、组织样本等），值是基因在样本中的表达量（需经过标准化，如量化标准化 QN）。

# 安装e1071包（如果需要用于支持向量机等算法）
# install.packages('e1071')

# 安装并加载BiocManager包，如果需要安装preprocessCore库
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("preprocessCore")

rm(list = ls())  # 清除工作空间
pacman::p_load(
  here,       # 用于处理文件路径
  data.table,  # 用于数据处理
  preprocessCore,  # 用于数据标准化
  e1071,  # 用于支持向量机等算法
  readr  # 用于读取CSV文件
)

here()

# 读取CSV文件
# inputFile = here("1_data/pancreatic_GEO/merged_GSE119794_GSE171485.csv")
# data <- read_csv(inputFile)
#
# # 设置输出文件路径
# outputFile = here("1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt")
#
# # 将数据写入TXT文件
# write.table(data, file = outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
inputFile = here("1_data/pancreatic_GEO/merged_GSE119794_GSE171485.txt")

# 加载包含CIBERSORT算法的自定义R脚本
source("2_src/00_GEOimmune.CIBERSORT.R")  # 通过CIBERSORT进行免疫细胞组成分析

# 使用CIBERSORT算法进行免疫细胞分布分析
# "ref.txt"是参考基因表达矩阵，用于解卷积免疫细胞成分
# inputFile是输入的基因表达数据
# perm = 1000 表示置换次数，QN=TRUE 表示是否进行量化标准化（Quantile Normalization）
outTab = CIBERSORT(here("1_data/ref.txt"), inputFile, perm = 1000, QN = TRUE)

# 筛选P-value小于1的结果
outTab = outTab[outTab[,"P-value"] < 1,]  # 保留p值小于1的有效样本

# 将数据转化为矩阵形式，去除最后三列无关数据
outTab = as.matrix(outTab[, 1:(ncol(outTab) - 3)])  # 去掉最后三列（通常是统计信息列）

# 添加列名并将结果写入文件
outTab = rbind(id = colnames(outTab), outTab)  # 在结果中添加列名
write.table(outTab, file = here("3_outputs/CIBERSORT-Results.txt"), sep = "\t", quote = FALSE, col.names = FALSE)  # 将结果保存为txt文件，当 quote = TRUE（默认值）时，函数会为字符型数据（如基因名、样本名）自动添加双引号，避免文本中包含空格、制表符等特殊字符时导致格式混乱。当 quote = FALSE 时，不会添加任何引号，输出的文本内容直接以原始形式保存（更符合常规的表格文件格式，便于后续分析工具读取）。若原始数据为 GeneA，quote = TRUE 会保存为 "GeneA"，而 quote = FALSE 会直接保存为 GeneA。col.names = FALSE 时，输出文件不会包含列名，直接从第一行开始写入数据内容。


