# 清理工作环境
rm(list = ls())

# 安装并加载所需包
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  here, TwoSampleMR, dplyr, tidyr
)

# --- 1. 初始化设置 ---

# 读取本地结局数据（替换原结局ID提取方式）
outcome_dat <- readRDS(here("3_outputs", "finngen_outcome_dat.rds"))
# 检查结局数据格式（确保包含TwoSampleMR所需列）
cat("本地结局数据包含", nrow(outcome_dat), "个SNP，关键列：",
    paste(colnames(outcome_dat), collapse = ", "), "\n")

# 读取暴露因子数据文件，获取唯一的暴露ID
exposure_data <- read.csv(here("1_data/imm.exposure.csv"))
exposureIDs <- unique(exposure_data$id.exposure)
cat("共检测", length(exposureIDs), "个暴露因子\n")

# 创建以今天日期命名的输出文件夹
today_date <- format(Sys.Date(), "%Y%m%d")
output_dir <- here("3_outputs", paste0(today_date, "_MR_with_local_outcome"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
print(paste("所有输出文件将保存到目录：", output_dir))

# 创建一个空的数据框，用于存储所有循环的最终结果
finalResults <- data.frame()


# --- 2. 循环执行孟德尔随机化分析 ---

# 循环处理每个暴露因子
for (i in seq_along(exposureIDs)) {
  exposureID <- exposureIDs[i]
  cat("------------------------------------------------------------------\n")
  cat("正在分析:", exposureID, "与本地结局数据的因果关系 (", i, "/", length(exposureIDs), ")\n")

  tryCatch({
    # A. 数据提取与协调

    # 提取暴露因子的工具变量 (p < 5e-8, 进行clump)
    exposure_dat <- extract_instruments(exposureID, p1 = 5e-8, clump = TRUE)
    if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
      stop("暴露因子没有满足条件的工具变量 (p < 5e-8)。")
    }
    cat("暴露因子", exposureID, "提取到", nrow(exposure_dat), "个工具变量\n")

    # 使用本地结局数据，无需从数据库提取
    # 直接基于暴露的SNP与本地结局数据匹配
    # （协调时会自动筛选匹配的SNP）
    dat <- harmonise_data(exposure_dat, outcome_dat)
    if (is.null(dat) || nrow(dat) == 0) {
      stop("暴露与结局数据协调后没有匹配的SNP。")
    }
    cat("协调后保留", nrow(dat), "个共有的SNP\n")

    # 保存协调后且通过筛选的SNP信息
    outTab <- dat[dat$mr_keep == TRUE, ]
    write.csv(outTab, file = file.path(output_dir, paste0("table.SNP_", exposureID, ".csv")), row.names = FALSE)

    # B. MR分析

    # 根据SNP数量选择合适的MR分析方法
    if (nrow(dat) == 1) {
      # 单个SNP时，仅使用Wald ratio方法
      mrResult <- mr(dat, method_list = "mr_wald_ratio")
    } else {
      # 多个SNP时，使用多种主流方法
      mrResult <- mr(dat, method_list = c("mr_ivw_mre", "mr_ivw_fe", "mr_egger_regression",
                                          "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
    }

    # 生成比值比（Odds Ratios）并保存结果
    mrTab <- generate_odds_ratios(mrResult)
    write.csv(mrTab, file = file.path(output_dir, paste0("table.MRresult_", exposureID, ".csv")), row.names = FALSE)

    # C. 敏感性分析 (仅当SNP数量 > 1时进行)

    if (nrow(dat) > 1) {
      # 异质性检验
      heterTab <- mr_heterogeneity(dat)
      write.csv(heterTab, file = file.path(output_dir, paste0("table.heterogeneity_", exposureID, ".csv")), row.names = FALSE)

      # 水平多效性检验
      pleioTab <- mr_pleiotropy_test(dat)
      write.csv(pleioTab, file = file.path(output_dir, paste0("table.pleiotropy_", exposureID, ".csv")), row.names = FALSE)
    } else {
      cat("SNP数量不足（仅", nrow(dat), "个），跳过异质性和多效性检验。\n")
    }

    # D. 结果汇总与可视化

    # 将本次循环的结果追加到最终结果数据框中
    finalResults <- rbind(finalResults, cbind(exposureID = exposureID, mrTab))

    # --- 绘图 ---
    # 对每个SNP单独进行MR分析，因果效应 = （SNP 与结局的效应 beta） / （SNP 与暴露的效应 beta）
    res_single <- mr_singlesnp(dat)

    # 绘制散点图
    pdf(file = file.path(output_dir, paste0("plot_", exposureID, "_1_scatter.pdf")), width = 7, height = 7)
    p1 <- mr_scatter_plot(mrResult, dat)
    print(p1[[1]])  # 提取图形对象
    dev.off()

    # 仅在多个SNP时绘制森林图、漏斗图和留一法图
    if (nrow(dat) > 1) {
      # 绘制森林图
      pdf(file = file.path(output_dir, paste0("plot_", exposureID, "_2_forest.pdf")), width = 7, height = 6)
      p2 <- mr_forest_plot(res_single)
      print(p2)
      dev.off()

      # 绘制漏斗图
      pdf(file = file.path(output_dir, paste0("plot_", exposureID, "_3_funnel.pdf")), width = 7, height = 7)
      p3 <- mr_funnel_plot(singlesnp_results = res_single)
      print(p3)
      dev.off()

      # 绘制留一法敏感性分析图
      pdf(file = file.path(output_dir, paste0("plot_", exposureID, "_4_leaveoneout.pdf")), width = 7, height = 6)
      p4 <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
      print(p4)
      dev.off()
    }

    cat("成功分析:", exposureID, "\n")

  }, error = function(e) {
    # 捕获错误并继续下一个循环
    cat("分析失败:", exposureID, "\n")
    cat("错误信息:", conditionMessage(e), "\n")
  })
}

# --- 3. 保存最终汇总结果 ---

write.csv(finalResults, file = file.path(output_dir, "summary_results_all_exposures.csv"), row.names = FALSE)

cat("------------------------------------------------------------------\n")
cat("全部分析完成！汇总结果已保存至：", output_dir, "\n")
