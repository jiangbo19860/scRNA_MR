# 安装和加载所需包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, grid, forestploter, pdftools)

rm(list = ls())  # 清除工作空间

# 确认项目根目录
here()

# 设置需展示的MR方法
selectMethod <- "Inverse variance weighted (fixed effects)"

# 定义目标目录并检查存在性
target_dir <- here("08_MRresult")
if (!dir.exists(target_dir)) {
  stop(paste("目录不存在:", target_dir))
}

# 读取并合并数据
files <- list.files(target_dir, pattern = "\\.csv$", full.names = TRUE)
print(paste("找到", length(files), "个CSV文件："))
head(files)

data <- data.frame()
for (i in files) {
  rt <- read.csv(i, header = TRUE, sep = ",", check.names = FALSE)
  data <- rbind(data, rt)
}
data <- data[data$method %in% selectMethod, ]
if (nrow(data) == 0) stop("过滤后的数据为空")

# --------------------------数据整理--------------------------
# 处理p值
data$pval_num <- as.numeric(data$pval)
data$pval <- ifelse(
  is.na(data$pval_num), "NA",
  ifelse(data$pval_num < 0.001, "<0.001", sprintf("%.3f", data$pval_num))
)

# 处理暴露变量（保留首个名称，其余空）
data$exposure_clean <- data$exposure
data[duplicated(data$exposure), "exposure_clean"] <- ""

# 处理SNP数量和空白列
data$nsnp <- as.character(data$nsnp)
data$nsnp[is.na(data$nsnp)] <- ""
data$space <- "                   "  # 空白列（放置可信区间）

# 格式化OR和95%CI
data$OR_95CI <- ifelse(
  is.na(data$or), "",
  sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95)
)

# --------------------------按100行拆分数据（每组含列名）--------------------------
# 1. 定义表头（列名和列数必须与数据子集完全一致）
header <- data.frame(
  exposure_clean = "Exposure",
  nsnp = "SNPs",
  method = "Method",
  pval = "P-value",
  space = "               ",
  OR_95CI = "OR (95% CI)",
  stringsAsFactors = FALSE
)

total_rows <- nrow(data)
rows_per_group <- 100  # 每100行一组
group_count <- ceiling(total_rows / rows_per_group)
group_indices <- list()

# 2. 生成每组数据行索引
for (g in 1:group_count) {
  start_row <- (g - 1) * rows_per_group + 1
  end_row <- min(g * rows_per_group, total_rows)
  group_indices[[g]] <- start_row:end_row
}

# --------------------------定义图形主题--------------------------
# 使用 forest_theme() 函数自定义森林图的外观主题
tm <- forest_theme(
  base_size = 14,  # 图形基础字体大小为14
  ci_pch = 16,     # 置信区间（CI）中点的形状，16表示实心圆
  ci_lty = 1,      # 置信区间线条的类型，1表示实线
  ci_lwd = 1.5,    # 置信区间线条的宽度为1.5
  ci_col = "black",# 置信区间线条的颜色为黑色
  ci_Theight = 0.2,# 置信区间线条末端"T"形的高度（误差线帽的大小）
  refline_gp = gpar(lty = "dashed", lwd = 1, col = "grey20"),  # 参考线（通常是无效线，如HR=1）的样式
  # 其中 gpar() 用于设置图形参数：lty为虚线，lwd为1，颜色为深灰色（grey20）
  xaxis_gp = gpar(cex = 0.9),  # x轴刻度文本的样式，cex=0.9表示比基础字体小10%
  footnote_gp = gpar(cex = 0.7, col = "blue")  # 脚注文本的样式，cex=0.7表示较小字体，颜色为蓝色
)

# --------------------------循环生成每组PDF（带列名）--------------------------
for (g in 1:group_count) {
  # 3. 提取当前组数据，并严格筛选与header一致的6列（关键修正）
  data_sub <- data[group_indices[[g]], c("exposure_clean", "nsnp", "method", "pval", "space", "OR_95CI")]

  # 4. 合并表头和数据（此时列数完全一致，可正常rbind）
  current_data <- rbind(header, data_sub)

  # 5. 构建绘图数据（确保列名与header一致）
  plot_data <- current_data  # 直接使用合并后的数据（已包含表头）
  rownames(plot_data) <- NULL  # 重置行名

  # 6. 准备效应值（表头行无效应值，用NA填充）
  est_values <- c(NA, data[group_indices[[g]], "or"])
  lower_values <- c(NA, data[group_indices[[g]], "or_lci95"])
  upper_values <- c(NA, data[group_indices[[g]], "or_uci95"])

  # 7. 计算X轴的动态范围
  # 取OR值95%置信区间下限的最小值，乘以0.95作为X轴左边界
  x_min <- min(data$or_lci95) * 0.95
  # 取OR值95%置信区间上限的最大值，乘以1.05作为X轴右边界
  x_max <- max(data$or_uci95) * 1.05

  # 8. 绘制森林图
  # 使用forest()函数创建森林图，将结果存储在plot对象中
  plot <- forest(
    data = plot_data,         # 指定用于绘图的数据集
    est = est_values,         # 指定效应值（如OR、HR等）所在的向量或列名
    lower = lower_values,     # 指定置信区间下限值所在的向量或列名
    upper = upper_values,     # 指定置信区间上限值所在的向量或列名
    ci_column = 5,            # 设置置信区间在图中显示的列位置（第5列）
    ref_line = 1,             # 添加参考线（通常用于森林图的无效线，如OR=1、HR=1）
    xlim = c(x_min, x_max),   # 使用动态计算的范围设置X轴
    theme = tm,               # 应用之前定义的自定义主题tm（控制整体外观）
    gp = gpar(cex = 1.1)      # 通过cex参数放大置信区间中点（1.5倍，可调整）
  )

  # 9. 美化图形（对置信区间样式进行调整，跳过表头行）
  boxcolor <- "red"  # 定义置信区间元素的填充颜色为红色

  # 循环遍历数据行：从第2行开始到最后一行（nrow(plot_data)获取总行数）
  # 跳过第1行是因为第1行通常是表头信息，不需要美化数据样式
  for (i in 2:nrow(plot_data)) {

    # 使用edit_plot()函数修改森林图中特定元素的样式
    plot <- edit_plot(
      plot,               # 指定要修改的森林图对象（即上面创建的plot）
      col = 5,            # 指定要修改的列（与ci_column保持一致，第5列的置信区间）
      row = i,            # 指定当前循环到的行（逐行修改数据行）
      which = "ci",       # 指定修改的元素类型为置信区间（ci = confidence interval）
      gp = gpar(fill = boxcolor, fontsize = 12)  # 设置图形参数：fill = boxcolor：用红色填充置信区间元素。fontsize = 12：设置该区域字体大小为12。
    )
  }

  # 10. p值<0.01的行加粗（跳过表头）
  # 注意：current_data的pval列第1行是表头"P-value"，需从第2行开始判断
  sig_logic <- c(FALSE, data[group_indices[[g]], "pval_num"] < 0.01 & !is.na(data[group_indices[[g]], "pval_num"]))
  sig_rows <- which(sig_logic)  # 显著行的索引（包含表头后的行）
  if (length(sig_rows) > 0) {
    for (i in sig_rows) {
      plot <- edit_plot(
        plot, col = 4, row = i,
        which = "text",
        gp = gpar(fontface = "bold")
      )
    }
  }

  # 11. 添加边框和居中对齐
  plot <- add_border(plot, part = "header", row = 1, where = "top", gp = gpar(lwd = 2))
  plot <- edit_plot(
    plot, col = 1:ncol(plot_data),
    which = "text",
    hjust = unit(0.5, "npc"),
    x = unit(0.5, "npc")
  )

  # 12. 输出PDF
  # 将文件保存到"results"子文件夹（需确保该文件夹已存在）
  pdf_file <- here("3_outputs", paste0("forest_group_", g, "_p_0.01.pdf"))
  pdf(pdf_file, width = 20, height = min(30, nrow(plot_data) * 0.8))
  print(plot)
  dev.off()

  # 13. 验证文件生成
  if (file.exists(pdf_file) && file.size(pdf_file) > 0) {
    message("已生成带列名的PDF：", pdf_file)
  } else {
    warning("PDF文件", pdf_file, "生成失败")
  }
}

# --------------------------提取p<0.01的结果并保存--------------------------
# 提取p<0.01的显著结果并保存到3_outputs文件夹
significant_data <- data[data$pval_num < 0.01 & !is.na(data$pval_num), ]
if (nrow(significant_data) > 0) {
  # 保存为CSV文件（使用here函数指定路径）
  write.csv(
    significant_data,
    here("3_outputs", "significant_results_p0.01.csv"),  # 保存到3_outputs文件夹
    row.names = FALSE,
    quote = FALSE
  )

  # 保存为TXT文件（使用here函数指定路径）
  write.table(
    significant_data[, c("exposure", "nsnp", "or", "or_lci95", "or_uci95", "pval")],
    here("3_outputs", "significant_results_p0.01.txt"),  # 保存到3_outputs文件夹
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  message("已保存显著结果（p<0.01）到3_outputs文件夹：", nrow(significant_data), "行")
} else {
  message("未发现p<0.01的显著结果")
}
