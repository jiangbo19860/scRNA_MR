# 输入：17步中得到的强相关的exposure_data_filtered.csv.

#devtools::install_github("mrcieu/gwasglue",force = TRUE)
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

setwd("/Users/lijiangbo/scRNA_MR/4_references/Stroke_MR+机器学习/37去除连锁不平衡性")     #设置工作目录

#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

#
exposure_data=read_exposure_data(filename="exposure_data_filtered.csv",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 pval_col = "pval.exposure",
                                 effect_allele_col="effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 samplesize_col = "samplesize.exposure",
                                 chr_col="chr.exposure", pos_col = "pos.exposure",
                                 clump=FALSE)
options(ieugwasr_api = 'https://gwas-api.mrcieu.ac.uk/')
# 注意去掉末尾的斜杠，保持简洁
#去除连锁不平衡的SNP，目的是使用 clump_data 函数对暴露数据进行 clumping 操作，以去除连锁不平衡的 SNP。
# 对暴露数据进行聚类，去除LD相关的SNP
exposure_dat_clumped = clump_data(
  exposure_data,  # 输入的暴露数据
  clump_kb = 10000,  # 聚类窗口：10000kb（10Mb），即仅考虑10Mb范围内的SNP间LD
  clump_r2 = 0.001   # LD阈值：当SNP间相关系数r² > 0.001时，视为存在LD并去除次要SNP
)
# 备选参数（更严格的阈值）：clump_kb=500, clump_r2=0.01
#exposure_dat_clumped=clump_data(exposure_data, clump_kb=500, clump_r2=0.01)
write.csv(exposure_dat_clumped, file="exposure_clumped.csv", row.names=F)

