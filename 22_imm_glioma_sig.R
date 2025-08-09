# MR结果forest与显著sig结果(固定效应逆方差加权法)，输出forest plot和sig_MR结果
rm(list = ls())
pacman::p_load(here, tidyverse, forestploter, grid)

# --------------------------创建输出文件夹（当天日期+任务名称）--------------------------
# 获取当天日期（格式：YYYYMMDD）
today_date <- format(Sys.Date(), "%Y%m%d")

# 定义输出文件夹名称（日期+核心任务）
output_dir_name <- paste0(today_date, "_MR_forest")

# 定义完整输出路径（3_outputs文件夹下）
base_output_dir <- here("3_outputs")
output_dir <- here(base_output_dir, output_dir_name)

# 创建文件夹（若已存在则不重复创建，避免覆盖）
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("所有输出将保存到：", output_dir)

# --------------------------读取数据--------------------------
# 设置原始MR结果文件所在路径
folder_path <- here("3_outputs/20250809_MR_with_local_outcome")

# 获取所有以table.MRresult_开头的CSV文件
files <- list.files(
  path = folder_path,
  pattern = "^table\\.MRresult_",
  full.names = TRUE
)

# 合并所有文件数据
data <- data.frame()
for (i in files) {
  rt <- read.csv(i, header = TRUE, sep = ",", check.names = FALSE)
  data <- rbind(data, rt)
}

# 筛选指定的MR方法（固定效应逆方差加权法）
selectMethod <- "Inverse variance weighted (fixed effects)"
data <- data[data$method %in% selectMethod, ]
if (nrow(data) == 0) stop("过滤后的数据为空，请检查MR方法名称是否正确")

# --------------------------数据整理--------------------------
# 处理p值（格式化显示）
data$pval_num <- as.numeric(data$pval)
data$pval <- ifelse(
  is.na(data$pval_num), "NA",
  ifelse(data$pval_num < 0.001, "<0.001", sprintf("%.3f", data$pval_num))
)

# 处理暴露变量（去重显示，仅保留首个名称）
data$exposure_clean <- data$exposure
data[duplicated(data$exposure), "exposure_clean"] <- ""

# 处理SNP数量和空白列（用于放置置信区间图形）
data$nsnp <- as.character(data$nsnp)
data$nsnp[is.na(data$nsnp)] <- ""
data$space <- "                   "

# 格式化OR值和95%置信区间（如1.234 (1.012 to 1.500)）
data$OR_95CI <- ifelse(
  is.na(data$or), "",
  sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95)
)

# --------------------------按100行拆分数据（避免单张图过大）--------------------------
# 定义表头（与数据列严格对应）
header <- data.frame(
  exposure_clean = "Exposure",
  nsnp = "SNPs",
  method = "Method",
  pval = "P-value",
  space = "               ",
  OR_95CI = "OR (95% CI)",
  stringsAsFactors = FALSE
)

# 计算分组数量和每组行索引
total_rows <- nrow(data)
rows_per_group <- 100
group_count <- ceiling(total_rows / rows_per_group)
group_indices <- list()

for (g in 1:group_count) {
  start_row <- (g - 1) * rows_per_group + 1
  end_row <- min(g * rows_per_group, total_rows)
  group_indices[[g]] <- start_row:end_row
}

# --------------------------定义森林图主题（外观样式）--------------------------
tm <- forest_theme(
  base_size = 14,
  ci_pch = 16,         # 置信区间中点形状（实心圆）
  ci_lty = 1,          # 置信区间线条类型（实线）
  ci_lwd = 1.5,        # 置信区间线条宽度
  ci_col = "black",    # 置信区间线条颜色
  ci_Theight = 0.2,    # 置信区间末端T形高度
  refline_gp = gpar(lty = "dashed", lwd = 1, col = "grey20"),  # 参考线（OR=1）样式
  xaxis_gp = gpar(cex = 0.9),  # x轴刻度字体大小
  footnote_gp = gpar(cex = 0.7, col = "blue"),  # 脚注样式
  mar = unit(c(5, 10, 10, 15), "mm")  # 上=5, 右=10, 下=15, 左=5（下边距增大）
)

# --------------------------循环生成森林图PDF（保存到目标文件夹）--------------------------
for (g in 1:group_count) {
  # 提取当前组数据并合并表头
  data_sub <- data[group_indices[[g]], c("exposure_clean", "nsnp", "method", "pval", "space", "OR_95CI")]
  current_data <- rbind(header, data_sub)
  plot_data <- current_data
  rownames(plot_data) <- NULL

  # 准备效应值（OR）及置信区间
  est_values <- c(NA, data[group_indices[[g]], "or"])
  lower_values <- c(NA, data[group_indices[[g]], "or_lci95"])
  upper_values <- c(NA, data[group_indices[[g]], "or_uci95"])

  # 计算X轴动态范围（基于置信区间）
  x_min <- min(data$or_lci95, na.rm = TRUE) * 0.95
  x_max <- max(data$or_uci95, na.rm = TRUE) * 1.05

  # 调用 forest 函数绘制森林图，将绘制结果赋值给 plot 对象
  plot <- forest(
    data = plot_data,  # 指定用于绘制森林图的数据集，包含表头及具体数据行，
    # 列应包含后续 est、lower、upper 等参数对应数据相关的列
    est = est_values,  # 传入效应值（如 OR 值等）的向量，表头行对应位置为 NA，数据行对应实际效应值，
    # 用于在森林图中展示每个条目的效应大小
    lower = lower_values,  # 传入效应值 95% 置信区间的下限向量，表头行对应位置为 NA，数据行是实际下限值，
    # 与 est、upper 配合绘制置信区间
    upper = upper_values,  # 传入效应值 95% 置信区间的上限向量，表头行对应位置为 NA，数据行是实际上限值
    ci_column = 5,        # 指定在数据集中用于展示置信区间图形的列索引（从 1 开始计数），
    # 即第 5 列会绘制置信区间的线条和点等图形元素
    ref_line = 1,         # 添加一条参考线，这里设置为 1，通常用于表示无效值（如 OR = 1 时代表无关联），
    # 帮助直观判断效应值与无效值的关系
    xlim = c(x_min, x_max),  # 设置 X 轴的范围，x_min 是根据数据计算的置信区间下限最小值再调整得到，
    # x_max 是置信区间上限最大值再调整得到，让图形在合适范围内展示
    theme = tm,  # 应用之前通过 forest_theme 函数定义好的森林图外观主题，
    # 控制字体、线条样式、颜色等整体视觉风格
    gp = gpar(cex = 1.1)  # 使用 gpar 函数设置图形参数，这里通过 cex 参数将图形中相关元素（如点、文本等）的大小放大为原来的 1.1 倍，让图形内容更清晰易看
  )

  # 美化置信区间（红色填充）
  boxcolor <- "red"
  for (i in 2:nrow(plot_data)) {  # 跳过表头行
    plot <- edit_plot(
      plot, col = 5, row = i, which = "ci",
      gp = gpar(fill = boxcolor, fontsize = 12)
    )
  }

  # p值<0.05的行加粗显示
  sig_logic <- c(FALSE, data[group_indices[[g]], "pval_num"] < 0.05 & !is.na(data[group_indices[[g]], "pval_num"]))
  sig_rows <- which(sig_logic)
  if (length(sig_rows) > 0) {
    for (i in sig_rows) {
      plot <- edit_plot(
        plot, col = 4, row = i, which = "text",
        gp = gpar(fontface = "bold")
      )
    }
  }

  # 添加表头边框和文本居中对齐
  plot <- add_border(plot, part = "header", row = 1, where = "top", gp = gpar(lwd = 2))
  plot <- edit_plot(
    plot, col = 1:ncol(plot_data), which = "text",
    hjust = unit(0.5, "npc"), x = unit(0.5, "npc")
  )

  # 保存森林图到目标文件夹
  pdf_file <- here(output_dir, paste0("forest_group_", g, "_p_0.05.pdf"))
  pdf(pdf_file, width = 20, height = min(30, nrow(plot_data) * 0.8))
  print(plot)
  dev.off()

  # 验证文件生成
  if (file.exists(pdf_file) && file.size(pdf_file) > 0) {
    message("已生成PDF：", pdf_file)
  } else {
    warning("PDF生成失败：", pdf_file)
  }
}

# --------------------------筛选、提取并保存显著结果（p<0.05）--------------------------
significant_data <- data[data$pval_num < 0.05 & !is.na(data$pval_num), ]
if (nrow(significant_data) > 0) {
  # 保存为CSV
  write.csv(
    significant_data,
    here(output_dir, "significant_results_p0.05.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  # 保存为TXT（制表符分隔）
  write.table(
    significant_data[, c("exposure", "nsnp", "or", "or_lci95", "or_uci95", "pval")],
    here(output_dir, "significant_results_p0.05.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  message("已保存显著结果到：", output_dir, "，共", nrow(significant_data), "行")
} else {
  message("未发现p<0.05的显著结果")
}
