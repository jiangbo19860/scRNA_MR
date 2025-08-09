# 用CIBERSORT算法分析GEO数据集中的免疫细胞的组成。

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
# "ref.txt" 是参考基因表达矩阵，用于解卷积免疫细胞成分
# inputFile 是输入的基因表达数据
# perm = 1000 表示置换次数，QN=TRUE 表示是否进行量化标准化（Quantile Normalization）
outTab = CIBERSORT(here("1_data/ref.txt"), inputFile, perm = 1000, QN = TRUE)

# 筛选P-value小于1的结果
outTab = outTab[outTab[,"P-value"] < 1,]  # 保留p值小于1的有效样本

# 将数据转化为矩阵形式，去除最后三列无关数据
outTab = as.matrix(outTab[, 1:(ncol(outTab) - 3)])  # 去掉最后三列（通常是统计信息列）

# 添加列名并将结果写入文件
outTab = rbind(id = colnames(outTab), outTab)  # 在结果中添加列名
write.table(outTab, file = here("3_outputs/CIBERSORT-Results.txt"), sep = "\t", quote = FALSE, col.names = FALSE)  # 将结果保存为txt文件
