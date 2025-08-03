rm(list = ls())

# 加载所需包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,
  tidyverse,
  data.table,
  dplyr,
  readr,
  ggplot2,
  TwoSampleMR
)

# 确认here()根目录（应返回/Users/lijiangbo/scRNA_MR）
cat("当前工作目录根路径：", here(), "\n\n")

# 读取包含id.exposure的CSV文件，提取唯一exposure ID
exposure_csv_path <- here("3_outputs", "significant_results_p0.05.csv")
if (!file.exists(exposure_csv_path)) {
  stop(paste0("文件 ", exposure_csv_path, " 不存在，请检查路径！"))
}
exposure_df <- read_csv(exposure_csv_path)
unique_exposures <- unique(exposure_df$id.exposure)
cat("提取到", length(unique_exposures), "个唯一的暴露ID：\n", paste(unique_exposures, collapse = "\n "), "\n\n")

# 定义结局ID
outcome_id <- "ebi-a-GCST90018893"

# 循环处理每个暴露ID
for (exposure_id in unique_exposures) {

  # 创建输出文件夹（暴露ID_结局ID）
  output_folder <- paste0(exposure_id, "_", outcome_id)
  output_dir <- here("3_outputs", output_folder)

  # 新建输出目录（若不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("已创建输出目录：", output_dir, "\n\n")
  } else {
    cat("输出目录已存在：", output_dir, "\n\n")
  }

  # 提取暴露数据
  bmi_exp_dat <- extract_instruments(outcomes = exposure_id)
  # 保存暴露工具变量
  write.csv(
    bmi_exp_dat,
    file = here(output_dir, "exposure_instruments.csv"),
    row.names = FALSE
  )
  cat("提取到", nrow(bmi_exp_dat), "个暴露工具变量SNP（暴露ID：", exposure_id, "）\n")

  # 提取结局数据
  chd_out_dat <- extract_outcome_data(
    snps = bmi_exp_dat$SNP,
    outcomes = outcome_id
  )
  # 保存结局数据
  write.csv(
    chd_out_dat,
    file = here(output_dir, "outcome_data.csv"),
    row.names = FALSE
  )
  cat("结局数据中匹配到", nrow(chd_out_dat), "个SNP（暴露ID：", exposure_id, "）\n")

  # 协调暴露和结局数据
  dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
  # 保存协调后的数据
  write.csv(
    dat,
    file = here(output_dir, "harmonised_data.csv"),
    row.names = FALSE
  )
  cat("数据协调后保留", nrow(dat), "个有效SNP（暴露ID：", exposure_id, "）\n")

  # 执行MR分析
  res <- mr(dat)
  # 保存MR分析结果
  write.csv(
    res,
    file = here(output_dir, "mr_results.csv"),
    row.names = FALSE
  )

  # 保存可用的MR方法列表
  mr_method_list() %>%
    write.csv(
      file = here(output_dir, "available_mr_methods.csv"),
      row.names = FALSE
    )

  # 敏感性分析：异质性检验
  heterogeneity_all <- mr_heterogeneity(dat)
  write.csv(
    heterogeneity_all,
    file = here(output_dir, "heterogeneity_test_all_methods.csv"),
    row.names = FALSE
  )

  heterogeneity_subset <- mr_heterogeneity(
    dat,
    method_list = c("mr_egger_regression", "mr_ivw")
  )
  write.csv(
    heterogeneity_subset,
    file = here(output_dir, "heterogeneity_test_egger_ivw.csv"),
    row.names = FALSE
  )

  # 水平多效性检验（Egger截距检验）
  pleiotropy_test <- mr_pleiotropy_test(dat)
  write.csv(
    pleiotropy_test,
    file = here(output_dir, "pleiotropy_test.csv"),
    row.names = FALSE
  )

  # 单SNP分析
  res_single <- mr_singlesnp(dat)
  res_single_fixed <- mr_singlesnp(dat, single_method = "mr_meta_fixed")
  res_single_ml <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")

  # 保存单SNP分析结果
  write.csv(
    res_single,
    file = here(output_dir, "single_snp_results_default.csv"),
    row.names = FALSE
  )
  write.csv(
    res_single_fixed,
    file = here(output_dir, "single_snp_results_fixed_effect.csv"),
    row.names = FALSE
  )
  write.csv(
    res_single_ml,
    file = here(output_dir, "single_snp_results_max_likelihood.csv"),
    row.names = FALSE
  )

  # 留一法敏感性分析
  res_loo <- mr_leaveoneout(dat)
  # 保存留一法结果
  write.csv(
    res_loo[[1]],  # 提取结果数据框
    file = here(output_dir, "leave_one_out_results.csv"),
    row.names = FALSE
  )

  # 绘制并保存散点图
  p1 <- mr_scatter_plot(res, dat)

  # 自定义散点图：添加x=0和y=0的黑线，并调整坐标轴样式
  p1_with_lines <- p1[[1]] +
    theme(
      # 控制坐标轴线条（x 轴和 y 轴）的样式，设置 linewidth 让线条更细
      axis.line.x.bottom = element_line(linewidth = 0.5, colour = "black"),
      axis.line.y = element_line(linewidth = 0.5, colour = "black"),
      # 调整坐标轴刻度线条粗细
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),

      # 调整坐标轴标题字体大小
      axis.title.x = element_text(vjust = 0, margin = margin(t = 5), size = 10),
      axis.title.y = element_text(angle = 90, vjust = 1, margin = margin(r = 5), size = 10),

      # 调整坐标轴刻度标签字体大小
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),

      # 调整图表标题字体大小
      plot.title = element_text(size = 12),

      # 调整图例字体大小
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold")
    )

  # 保存修改后的散点图
  ggsave(
    plot = p1_with_lines,
    file = here(output_dir, "mr_scatter_plot_with_axes.pdf"),
    width = 7, height = 7
  )

  # 仅用Egger和IVW方法绘制散点图
  res_subset <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
  p1_subset <- mr_scatter_plot(res_subset, dat)
  ggsave(
    plot = p1_subset[[1]],
    file = here(output_dir, "mr_scatter_plot_egger_ivw.pdf"),
    width = 7, height = 7
  )

  # 绘制并保存森林图
  p2 <- mr_forest_plot(res_single)
  ggsave(
    plot = p2[[1]],
    file = here(output_dir, "forest_plot_default.pdf"),
    width = 10, height = 8 + 0.1 * nrow(dat)  # 高度随SNP数量调整
  )

  # 用指定方法绘制森林图
  res_single_subset <- mr_singlesnp(
    dat,
    all_method = c("mr_ivw", "mr_two_sample_ml")
  )
  p2_subset <- mr_forest_plot(res_single_subset)
  ggsave(
    plot = p2_subset[[1]],
    file = here(output_dir, "forest_plot_ivw_ml.pdf"),
    width = 10, height = 8 + 0.1 * nrow(dat)
  )

  # 绘制并保存留一法图
  p3 <- mr_leaveoneout_plot(res_loo)
  ggsave(
    plot = p3[[1]],
    file = here(output_dir, "leave_one_out_plot.pdf"),
    width = 10, height = 6
  )

  # 保存Egger方法的留一法结果（前几行）
  mr_leaveoneout(dat, method = mr_egger_regression) %>%
    head() %>%
    write.csv(
      file = here(output_dir, "leave_one_out_egger_head.csv"),
      row.names = FALSE
    )

  # 绘制并保存漏斗图
  p4 <- mr_funnel_plot(res_single)
  ggsave(
    plot = p4[[1]],
    file = here(output_dir, "funnel_plot.pdf"),
    width = 7, height = 7
  )

  cat("暴露ID", exposure_id, "分析完成！结果已保存至：", output_dir, "\n\n")
}

cat("所有暴露ID分析完成！\n")
