## https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html

rm(list = ls())
here()
pacman::p_load(
  here,
  tidyverse,
  data.table,
  dplyr,
  readr,
  ggplot2,
  TwoSampleMR
)
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')
dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
res <- mr(dat)
mr_method_list()

# Sensitivity analyses
mr_heterogeneity(dat)
mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
mr_pleiotropy_test(dat) # 水平多效性检验
res_single <- mr_singlesnp(dat)
res_single <- mr_singlesnp(dat, single_method = "mr_meta_fixed")
res_single <- mr_singlesnp(dat, all_method = "mr_two_sample_ml")
res_loo <- mr_leaveoneout(dat)
head(res_loo)
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
length(p1)
ggsave(p1[[1]], file = "filename.pdf", width = 7, height = 7)
# 仅绘制 MR Egger 和 IVW 因果效应估计的线条：
res <- mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
p2 <- mr_forest_plot(res_single)
p2[[1]]
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw", "mr_two_sample_ml"))
p2 <- mr_forest_plot(res_single)
p2[[1]]

p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
head(mr_leaveoneout(dat, method = mr_egger_regression))
p4 <- mr_funnel_plot(res_single)
p4[[1]]

