rm(list = ls())  # 清除工作空间
# 引用包
library(TwoSampleMR)

outcome <- read.csv("~/scRNA_MR/4_references/Stroke_MR+机器学习/09对有意义的暴露因素单独出图表/outcome.csv")
library(readxl)
myoutcome <- read_excel("~/scRNA_MR/3_outputs/有显著因果关系的25个暴露的SNPs合并后.xlsx")

colnames(outcome)
colnames(myoutcome)


exposureFile <- "imm.exposure.F.csv"  # 暴露数据文件
outcomeFile <- "outcome.csv"  # 结局数据文件
sigGutFile <- "immune-disease.IVWfilter.csv"  # IVW方法过滤的结果文件
outcomeName <- "Sepsis"  # 结局中展示疾病的名字
setwd("H:\\我的教程\\MR+机器学习\\09对有意义的暴露因素单独出图表")  # 设置工作目录

# 读取暴露数据输入文件
rt <- read.csv(exposureFile, header = TRUE, sep = ",", check.names = FALSE)

# 读取IVW方法过滤的结果文件
sigGut <- read.csv(sigGutFile, header = TRUE, sep = ",", check.names = FALSE)

# 对疾病相关的进行循环
for (i in sigGut$exposure) {
  # 替换特殊字符为下划线
  sanitized_name <- gsub("[^A-Za-z0-9_]", "_", i)
  singleExposureFile <- paste0(sanitized_name, ".exposure.csv")

  exposure_set <- rt[rt$exposure == i, ]
  write.csv(exposure_set, file = singleExposureFile, row.names = FALSE)

  # 读取暴露数据
  exposure_dat <- read_exposure_data(filename = singleExposureFile,
                                     sep = ",",
                                     snp_col = "SNP",
                                     beta_col = "beta.exposure",
                                     se_col = "se.exposure",
                                     pval_col = "pval.exposure",
                                     effect_allele_col = "effect_allele.exposure",
                                     other_allele_col = "other_allele.exposure",
                                     eaf_col = "eaf.exposure",
                                     phenotype_col = "exposure",
                                     samplesize_col = "samplesize.exposure",
                                     chr_col = "chr.exposure", pos_col = "pos.exposure",
                                     clump = FALSE)

  # 读取结局数据
  outcome_data <- read_outcome_data(snps = exposure_dat$SNP,
                                    filename = outcomeFile, sep = ",",
                                    snp_col = "SNP",
                                    beta_col = "beta.outcome",
                                    se_col = "se.outcome",
                                    effect_allele_col = "effect_allele.outcome",
                                    other_allele_col = "other_allele.outcome",
                                    pval_col = "pval.outcome",
                                    eaf_col = "eaf.outcome")

  # 将暴露数据和结局数据合并
  outcome_data$outcome <- outcomeName
  dat <- harmonise_data(exposure_dat, outcome_data)

  # 输出用于孟德尔随机化的工具变量
  outTab <- dat[dat$mr_keep == "TRUE", ]
  write.csv(outTab, file = paste0(sanitized_name, ".table.SNP.csv"), row.names = FALSE)

  # MR-PRESSO异常值检测(偏倚的SNP)
  presso <- run_mr_presso(dat)
  write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = paste0(sanitized_name, ".table.MR-PRESSO_Global.csv"))
  write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = paste0(sanitized_name, ".table.MR-PRESSO_Outlier.csv"))

  # 孟德尔随机化分析
  mrResult <- mr(dat)

  # 对结果进行OR值的计算
  mrTab <- generate_odds_ratios(mrResult)
  write.csv(mrTab, file = paste0(sanitized_name, ".table.MRresult.csv"), row.names = FALSE)

  # 异质性分析
  heterTab <- mr_heterogeneity(dat)
  write.csv(heterTab, file = paste0(sanitized_name, ".table.heterogeneity.csv"), row.names = FALSE)

  # 多效性检验
  pleioTab <- mr_pleiotropy_test(dat)
  write.csv(pleioTab, file = paste0(sanitized_name, ".table.pleiotropy.csv"), row.names = FALSE)

  # 绘制散点图
  pdf(file = paste0(sanitized_name, ".scatter_plot.pdf"), width = 7, height = 6.5)
  p1 <- mr_scatter_plot(mrResult, dat)
  print(p1)
  dev.off()

  # 森林图
  res_single <- mr_singlesnp(dat)  # 得到每个工具变量对结局的影响
  pdf(file = paste0(sanitized_name, ".forest.pdf"), width = 6.5, height = 5)
  p2 <- mr_forest_plot(res_single)
  print(p2)
  dev.off()

  # 漏斗图
  pdf(file = paste0(sanitized_name, ".funnel_plot.pdf"), width = 6.5, height = 6)
  p3 <- mr_funnel_plot(singlesnp_results = res_single)
  print(p3)
  dev.off()

  # 留一法敏感性分析
  pdf(file = paste0(sanitized_name, ".leaveoneout.pdf"), width = 6.5, height = 5)
  p4 <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
  print(p4)
  dev.off()
}
