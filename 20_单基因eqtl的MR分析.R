# 要先从'/Users/lijiangbo/scRNA_MR/1_data/gwas extracted_data 翻译.csv'中第4列查询基因名称, 然后到https://gwas.mrcieu.ac.uk中搜索并下载vcf文件.
# 输出：工具变量表格（table.SNP.csv）；MR 因果效应结果（table.MRresult.csv）；异质性、多效性检验结果（table.heterogeneity.csv、table.pleiotropy.csv）；散点图，森林图，漏斗图，留一法敏感性分析图。

#devtools::install_github("mrcieu/gwasglue",force = TRUE)
#BiocManager::install("VariantAnnotation")

#install.packages("htmltools")
#install.packages("devtools")
#library(devtools)
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
#usethis::edit_r_environ()
#GITHUB_TOKEN="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"


#1 ，GITHUB_TOKEN  2.网络（tizi）  3.不要更新包，跳过  4.R版本 4.3.0以上


rm(list = ls())
library(VariantAnnotation)  # 处理VCF格式的GWAS数据
library(gwasglue)           # 转换GWAS数据格式，适配MR分析
library(TwoSampleMR)        # 核心MR分析工具包

exposureFile="exposure.F.csv"            #暴露数据
outcomeFile="ieu-b-4980.vcf.gz" # 结局的遗传关联数据。
outcomeName="SPESIS"     # 图形中展示疾病的名称
setwd("/Users/lijiangbo/scRNA_MR/4_references/Stroke_MR+机器学习/39单基因与疾病的MR分析")     #设置工作目录

#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#读取结局数据的vcf文件,并对数据进行格式转换
vcfRT=readVcf(outcomeFile)
# 将VCF格式转换为MR分析所需的格式（包含SNP与结局的关联信息）
outcomeData=gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="outcome")

#从结局数据中提取工具变量
outcomeTab=merge(exposure_dat, outcomeData, by.x="SNP", by.y="SNP")
write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv")

#读取整理好的结局数据
outcome_data=read_outcome_data(snps=exposure_dat$SNP,
                               filename="outcome.csv", sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta.outcome",
                               se_col = "se.outcome",
                               effect_allele_col = "effect_allele.outcome",
                               other_allele_col = "other_allele.outcome",
                               pval_col = "pval.outcome",
                               eaf_col = "eaf.outcome")

#将暴露数据和结局数据合并
outcome_data$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcome_data)

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#孟德尔随机化分析
mrResult=mr(dat)

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="3_outputs/table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="3_outputs/table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="3_outputs/table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="3_outputs/pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="3_outputs/pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="3_outputs/pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="3_outputs/pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()

