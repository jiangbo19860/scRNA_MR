# 从 eQTL（表达数量性状基因座）数据的 VCF 文件中提取遗传变异（SNP）信息，并筛选出与目标基因表达显著相关的强关联 SNP.

#devtools::install_github("mrcieu/gwasglue",force = TRUE)
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")

rm(list = ls())  # 清空工作空间
#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
setwd("/Users/lijiangbo/scRNA_MR/4_references/Stroke_MR+机器学习/36筛选强关联性SNP")

exposureFile="eqtl-a-ENSG00000015133.vcf.gz"

#读取的vcf文件,并对数据进行格式转换
vcfRT=readVcf(exposureFile)
exposureData=gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
colnames(exposureData)

write.csv(exposureData, file="exposureFile.csv")

# 保留与基因表达显著关联的SNP（P值 < 5e-06，即5×10^-6，属于严格的显著性阈值）
exposure_data_filtered <- exposureData[exposureData$pval.exposure < 5e-06, ]

# 将筛选后的强关联SNP数据保存为CSV文件
write.csv(exposure_data_filtered, file = "exposure_data_filtered.csv", row.names = FALSE)


