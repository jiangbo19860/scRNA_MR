rm(list = ls())

here()

# 载入 TwoSampleMR 包
library(TwoSampleMR)

# 读取暴露因子的ID数据
exposureIDs <- read.table("extracted_data.txt", header = FALSE)$V1

outcomeID <- "ebi-a-GCST90018893"    # 结果数据ID

# 创建一个空的数据框，用于存储所有结果
finalResults <- data.frame()

# 循环处理每个暴露因子
for (i in seq_along(exposureIDs)) {
  exposureID <- exposureIDs[i]
  cat("Analysing", exposureID, "on", outcomeID, "(", i, "of", length(exposureIDs), ")\n")

  tryCatch({
    # 提取暴露因子数据
    exposure_dat <- extract_instruments(exposureID, p1 = 1e-8, clump = TRUE)

    if (nrow(exposure_dat) == 0) {
      stop("Not enough SNPs available for analysis")
    }

    # 提取结果数据
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcomeID)

    if (nrow(outcome_dat) == 0) {
      stop("Not enough SNPs available for Heterogeneity or pleiotropy analysis")
    }

    # 将暴露因子数据和结果数据进行协调
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # 提取需要保留的数据用于后续分析
    outTab <- dat[dat$mr_keep == "TRUE", ]
    write.csv(outTab, file = paste0("table.SNP_", exposureID, ".csv"), row.names = FALSE)

    # 进行 MR 分析
    mrResult <- mr(dat)

    # 选择特定的 MR 方法进行分析
    mrResult <- mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))

    # 生成 Odds Ratios（OR）的表格
    mrTab <- generate_odds_ratios(mrResult)
    write.csv(mrTab, file = paste0("table.MRresult_", exposureID, ".csv"), row.names = FALSE)

    # 进行异质性检验
    if (nrow(dat) > 1) {
      heterTab <- mr_heterogeneity(dat)
      write.csv(heterTab, file = paste0("table.heterogeneity_", exposureID, ".csv"), row.names = FALSE)
    } else {
      cat("Not enough SNPs available for Heterogeneity analysis of", exposureID, "on", outcomeID, "\n")
    }

    # 进行多重性检验
    pleioTab <- mr_pleiotropy_test(dat)
    if (!is.null(pleioTab)) {
      write.csv(pleioTab, file = paste0("table.pleiotropy_", exposureID, ".csv"), row.names = FALSE)
    } else {
      cat("Not enough SNPs available for pleiotropy analysis of", exposureID, "on", outcomeID, "\n")
    }

    # 将每次循环中生成的结果追加到总的结果数据框
    finalResults <- rbind(finalResults, cbind(exposureID = exposureID, mrTab))

    # 绘制散点图
    pdf(file = paste0(exposureID, "_scatter_plot.pdf"), width = 7, height = 6.5)
    p1 <- mr_scatter_plot(mrResult, dat)
    print(p1)
    dev.off()

    # 森林图
    res_single <- mr_singlesnp(dat)      # 得到每个工具变量对结局的影响
    pdf(file = paste0(exposureID, "_forest.pdf"), width = 6.5, height = 5)
    p2 <- mr_forest_plot(res_single)
    print(p2)
    dev.off()

    # 漏斗图
    pdf(file = paste0(exposureID, "_funnel_plot.pdf"), width = 6.5, height = 6)
    p3 <- mr_funnel_plot(singlesnp_results = res_single)
    print(p3)
    dev.off()

    # 留一法敏感性分析
    pdf(file = paste0(exposureID, "_leaveoneout.pdf"), width = 6.5, height = 5)
    p4 <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    print(p4)
    dev.off()

  }, error = function(e) {
    # 捕获错误并继续下一个暴露因子
    cat("Failed to analyse exposureID:", exposureID, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
}

# 生成汇总表格
write.csv(finalResults, file = "summary_results_Mendelian_randomization.csv", row.names = FALSE)
