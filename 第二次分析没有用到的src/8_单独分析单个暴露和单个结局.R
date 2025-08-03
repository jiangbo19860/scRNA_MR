rm(list = ls())

# 加载所需包（新增MR-PRESSO）
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,
  tidyverse,
  data.table,
  dplyr,
  readr,
  ggplot2,
  TwoSampleMR,
  MRPRESSO  # 加载MR-PRESSO包
)

# 确认here()根目录
cat("当前工作目录根路径：", here(), "\n\n")

# 定义暴露和结局ID
exposure_id <- "ebi-a-GCST90001670"
outcome_id <- "ebi-a-GCST90018893"

# 创建输出文件夹
output_folder <- paste0(exposure_id, "_", outcome_id)
output_dir <- here("3_outputs", output_folder)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("已创建输出目录：", output_dir, "\n\n")
} else {
  cat("输出目录已存在：", output_dir, "\n\n")
}

# 提取暴露数据
exp_dat <- extract_instruments(outcomes = exposure_id)
write.csv(
  exp_dat,
  file = here(output_dir, "exposure_instruments.csv"),
  row.names = FALSE
)
cat("提取到", nrow(exp_dat), "个暴露工具变量SNP\n")

# 提取结局数据
out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = outcome_id
)
write.csv(
  out_dat,
  file = here(output_dir, "outcome_data.csv"),
  row.names = FALSE
)
cat("结局数据中匹配到", nrow(out_dat), "个SNP\n")

# 协调暴露和结局数据
dat <- harmonise_data(exp_dat, out_dat)
dat_keep <- dat[dat$mr_keep == "TRUE", ]  # 保留有效SNP
snp_count <- nrow(dat_keep)
write.csv(
  dat_keep,
  file = here(output_dir, "harmonised_data.csv"),
  row.names = FALSE
)
cat("数据协调后保留", snp_count, "个有效SNP\n")

# 执行MR分析
res <- mr(dat_keep)
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
heterogeneity_all <- mr_heterogeneity(dat_keep)
write.csv(
  heterogeneity_all,
  file = here(output_dir, "heterogeneity_test_all_methods.csv"),
  row.names = FALSE
)

heterogeneity_subset <- mr_heterogeneity(
  dat_keep,
  method_list = c("mr_egger_regression", "mr_ivw")
)
write.csv(
  heterogeneity_subset,
  file = here(output_dir, "heterogeneity_test_egger_ivw.csv"),
  row.names = FALSE
)

# 水平多效性检验（Egger截距检验）
if (snp_count >= 3) {  # Egger需要至少3个SNP
  pleiotropy_test <- mr_pleiotropy_test(dat_keep)
  write.csv(
    pleiotropy_test,
    file = here(output_dir, "pleiotropy_test.csv"),
    row.names = FALSE
  )
} else {
  cat("SNP数量不足3个，跳过Egger多效性检验\n")
}

# 新增：MR-PRESSO分析（需至少3个SNP）
if (snp_count >= 3) {
  cat("执行MR-PRESSO分析...\n")
  tryCatch({
    # 运行MR-PRESSO（置换次数1000，显著性阈值0.05）
    presso_result <- run_mr_presso(
      dat = dat_keep,
      NbDistribution = 1000,
      SignifThreshold = 0.05
    )

    # 保存全局检验结果（评估整体多效性）
    write.csv(
      presso_result[[1]]$`MR-PRESSO results`$`Global Test`,
      file = here(output_dir, "mr_presso_global_test.csv"),
      row.names = FALSE
    )

    # 保存异常值检验结果（识别具体异常SNP）
    write.csv(
      presso_result[[1]]$`MR-PRESSO results`$`Outlier Test`,
      file = here(output_dir, "mr_presso_outlier_test.csv"),
      row.names = FALSE
    )

    # 保存校正异常值后的因果效应结果
    write.csv(
      presso_result[[1]]$`MR-PRESSO results`$`Causal Estimate`,
      file = here(output_dir, "mr_presso_causal_estimate.csv"),
      row.names = FALSE
    )

    cat("MR-PRESSO分析完成，结果已保存\n")
  }, error = function(e) {
    cat("MR-PRESSO分析失败：", e$message, "\n")
  })
} else {
  cat("SNP数量不足3个，跳过MR-PRESSO分析\n")
}

# 单SNP分析
res_single <- mr_singlesnp(dat_keep)
res_single_fixed <- mr_singlesnp(dat_keep, single_method = "mr_meta_fixed")
res_single_ml <- mr_singlesnp(dat_keep, all_method = "mr_two_sample_ml")

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
res_loo <- mr_leaveoneout(dat_keep)
write.csv(
  res_loo[[1]],
  file = here(output_dir, "leave_one_out_results.csv"),
  row.names = FALSE
)

# 绘制并保存散点图（自定义样式）
p1 <- mr_scatter_plot(res, dat_keep)
p1_with_lines <- p1[[1]] +
  theme(
    axis.line.x.bottom = element_line(linewidth = 0.5, colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    axis.title.x = element_text(vjust = 0, margin = margin(t = 5), size = 10),
    axis.title.y = element_text(angle = 90, vjust = 1, margin = margin(r = 5), size = 10),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )
ggsave(
  plot = p1_with_lines,
  file = here(output_dir, "mr_scatter_plot_with_axes.pdf"),
  width = 7, height = 7
)

# 仅用Egger和IVW方法绘制散点图（若SNP≥3）
if (snp_count >= 3) {
  res_subset <- mr(dat_keep, method_list = c("mr_egger_regression", "mr_ivw"))
  p1_subset <- mr_scatter_plot(res_subset, dat_keep)
  ggsave(
    plot = p1_subset[[1]],
    file = here(output_dir, "mr_scatter_plot_egger_ivw.pdf"),
    width = 7, height = 7
  )
} else {
  cat("SNP数量不足3个，跳过Egger+IVW散点图\n")
}

# 绘制并保存森林图
p2 <- mr_forest_plot(res_single)
ggsave(
  plot = p2[[1]],
  file = here(output_dir, "forest_plot_default.pdf"),
  width = 10, height = 8 + 0.1 * snp_count
)

# 用指定方法绘制森林图
res_single_subset <- mr_singlesnp(
  dat_keep,
  all_method = if (snp_count >= 3) c("mr_ivw", "mr_two_sample_ml") else "mr_ivw"
)
p2_subset <- mr_forest_plot(res_single_subset)
ggsave(
  plot = p2_subset[[1]],
  file = here(output_dir, "forest_plot_ivw_ml.pdf"),
  width = 10, height = 8 + 0.1 * snp_count
)

# 绘制并保存留一法图
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(
  plot = p3[[1]],
  file = here(output_dir, "leave_one_out_plot.pdf"),
  width = 10, height = 6
)

# 保存Egger方法的留一法结果（仅当SNP≥3）
if (snp_count >= 3) {
  mr_leaveoneout(dat_keep, method = mr_egger_regression) %>%
    head() %>%
    write.csv(
      file = here(output_dir, "leave_one_out_egger_head.csv"),
      row.names = FALSE
    )
}

# 绘制并保存漏斗图（需至少3个SNP）
if (snp_count >= 3) {
  p4 <- mr_funnel_plot(res_single)
  ggsave(
    plot = p4[[1]],
    file = here(output_dir, "funnel_plot.pdf"),
    width = 7, height = 7
  )
} else {
  cat("SNP数量不足3个，跳过漏斗图\n")
}

cat("所有分析完成！结果已保存至：", output_dir, "\n")
