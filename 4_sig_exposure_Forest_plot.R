rm(list = ls())  # 清空工作空间

# 安装和加载所需包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(grid, readr, dplyr, stringr, fs, forestploter, here)

# 1. 安装和加载所需包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(grid, readr, dplyr, stringr, fs, forestploter, here)

# 1. 设置路径和创建输出文件夹（核心修改：自动生成带序号的新文件夹）
today_date <- format(Sys.Date(), "%Y%m%d")
base_dir <- here("3_outputs")
base_name <- paste0(today_date, "_sig_MRresults")
output_dir <- here(base_dir, base_name)

# 检查文件夹是否存在，若存在则添加序号
counter <- 1
while (dir_exists(output_dir)) {
  output_dir <- here(base_dir, paste0(base_name, "_", counter))  # 如20250803_sig_MRresults_1
  counter <- counter + 1
}

# 创建最终的输出文件夹
dir_create(output_dir)
message("输出文件夹已创建: ", output_dir)

# 2. 从significant_results_p0.05.csv提取id.exposure的唯一值
sig_results_path <- here("3_outputs/20250802_MR_forest/significant_results_p0.05.csv")
sig_results <- read_csv(sig_results_path, show_col_types = FALSE)

# 提取并处理唯一的id.exposure值
target_ids <- sig_results %>%
  pull(id.exposure) %>%
  str_trim() %>%
  unique() %>%
  na.omit()

if (length(target_ids) == 0) {
  stop("在significant_results_p0.05.csv中未找到有效的id.exposure值")
}
message("提取到的id.exposure唯一值数量: ", length(target_ids))
message("部分id示例: ", paste(head(target_ids), collapse = ", "))

# 3. 从20250801_MR文件夹筛选符合条件的文件
source_dir <- here("3_outputs/20250801_MR")
all_csv_files <- dir_ls(path = source_dir, type = "file", regexp = "\\.csv$")

# 构建匹配模式：以table.MRresult_开头且包含目标id
pattern <- paste0("^table\\.MRresult_.*(", paste(target_ids, collapse = "|"), ")")
matched_files <- all_csv_files[str_detect(basename(all_csv_files), pattern)]

if (length(matched_files) == 0) {
  stop("未找到符合条件的文件（以table.MRresult_开头且包含目标id）")
}
message("找到符合条件的文件数量: ", length(matched_files))

# 复制文件到输出文件夹
file_copy(matched_files, output_dir, overwrite = TRUE)
message("已将符合条件的文件复制到输出文件夹")

# 4. 读取并合并数据（根据文件行数选择对应方法）
mr_files <- dir_ls(output_dir, type = "file", regexp = "\\.csv$")
data <- data.frame()

# 循环合并所有文件，根据文件行数选择方法
for (file in mr_files) {
  rt <- read_csv(file, show_col_types = FALSE)
  n_rows <- nrow(rt)  # 获取当前文件的行数

  # 根据行数筛选对应方法
  if (n_rows == 1) {
    # 单行数据：选择Wald ratio方法
    rt_filtered <- rt %>% filter(method == "Wald ratio")
    method_used <- "Wald ratio"
  } else if (n_rows >= 2) {
    # ≥2行数据：选择Inverse variance weighted (fixed effects)方法
    rt_filtered <- rt %>% filter(method == "Inverse variance weighted (fixed effects)")
    method_used <- "Inverse variance weighted (fixed effects)"
  } else {
    # 空文件跳过
    message("警告：文件", basename(file), "为空，已跳过")
    next
  }

  # 检查筛选后的数据是否存在
  if (nrow(rt_filtered) == 0) {
    message("警告：文件", basename(file), "中未找到", method_used, "方法的数据，已跳过")
    next
  }

  # 合并数据
  data <- rbind(data, rt_filtered)
}

# 检查合并后的数据是否为空
if (nrow(data) == 0) {
  stop("所有文件中均未找到符合条件的方法数据（Wald ratio或Inverse variance weighted (fixed effects)）")
}
message("合并后保留的数据行数: ", nrow(data))

# 5. 数据处理与格式化（确保添加空格列）
data <- data %>%
  # 先添加空格列（关键修正：确保该列存在）
  mutate(` ` = strrep(" ", 15)) %>%  # 15个空格，用于分隔表格与置信区间
  # 格式化OR值和95%置信区间
  mutate(
    `OR(95% CI)` = ifelse(
      is.na(or),
      "",
      sprintf("%.3f (%.3f to %.3f)", or, or_lci95, or_uci95)
    ),
    # 美化P值展示
    pval = ifelse(
      pval < 0.001,
      "<0.001",
      sprintf("%.3f", pval)
    ),
    # 处理暴露名称缺失值
    exposure = ifelse(is.na(exposure), "", exposure)
  ) %>%
  # 去除重复的暴露名称（仅保留首个出现的名称）
  mutate(exposure = ifelse(duplicated(exposure), "", exposure))

# 检查列是否存在（调试用）
print("数据框的列名：")
print(colnames(data))  # 确认包含 " " 列
new_headers <- c("Exposure", "NSNP", "Method", "P_val", " ", "OR (95% CI)")
plot_data <- data[, c("exposure", "nsnp", "method", "pval", " ", "OR(95% CI)")]
colnames(plot_data) <- new_headers
colnames(plot_data)
plot_data <- plot_data %>%
  mutate(
    Method = str_replace(
      Method,
      "Inverse variance weighted \\(fixed effects\\)",  # 匹配原始字符串（注意转义的转义）
      "Inverse variance weighted"  # 替换为简化后的字符串
    )
  )
write.csv(plot_data, file = here(output_dir, "processed_sig_exposure_forest.csv"), row.names = FALSE)

# 6. 准备森林图参数
lineVec <- cumsum(c(1, table(data$exposure[data$exposure != ""])))

# 定义森林图的外观主题（使用新参数替代已废弃的旧参数，确保代码兼容性）
tm <- forest_theme(  # forest_theme()是 forestploter 包中用于统一设置森林图外观的核心函数，
  base_size = 18,  # 全局基础字体大小（影响表头、正文等默认字体，单位为磅）
  ci_pch = 16,     # 置信区间中点符号的类型（16表示实心圆●，1=空心圆○，15=实心正方形■）
  ci_lty = 1,      # 置信区间线段的线型（1表示实线，2表示虚线，3表示点线等）
  ci_lwd = 1.5,    # 置信区间线段的线宽（数值越大线越粗，1.5为适中宽度）
  ci_col = "black",# 置信区间线段和中点的颜色（"black"表示黑色，可替换为"red"、"#FF0000"等）
  ci_Theight = 0.2,# 置信区间线段两端"T"形竖线的高度（0.2表示相对高度，数值越小竖线越短）
  refline_gp = gpar(  # 参考线（通常为OR=1的竖线）的样式设置
    lty = "dashed",  # 参考线线型（"dashed"表示虚线）
    lwd = 1,         # 参考线线宽（1为默认细实线）
    col = "grey20"   # 参考线颜色（"grey20"表示深灰色，避免与置信区间线冲突）
  ),
  xaxis_gp = gpar(cex = 0.8),  # X轴刻度文本的样式设置（cex=0.8表示相对基础字体缩小为80%）
  footnote_gp = gpar(          # 脚注文本的样式设置
    cex = 0.6,                 # 脚注字体相对大小（0.6表示基础字体的60%）
    col = "blue"               # 脚注颜色（"blue"表示蓝色，用于区分正文）
  )
)

# 7. 绘制森林图（使用处理好的plot_data作为输入数据，确保列名与新表头new_headers完全匹配，避免表头显示异常）
plot <- forestploter::forest(
  data = plot_data,  # 核心参数：指定用于绘图的数据集（已处理列名和内容的plot_data）
  est = data$or,     # 指定效应量（如OR值）的列，用于绘制置信区间的中点
  lower = data$or_lci95,  # 指定效应量95%置信区间的下限值列
  upper = data$or_uci95,  # 指定效应量95%置信区间的上限值列
  ci_column = 5,     # 指定在第5列绘制置信区间线段（对应数据集中的空格分隔列，视觉上更清晰）
  ref_line = 1,      # 添加参考线（通常为效应量=1的竖线，用于判断效应方向：>1为正向，<1为负向）
  xlim = c(0, max(data$or_uci95, na.rm = TRUE) + 0.5),  # 设置X轴范围：从0到最大置信区间上限+0.5（避免线段超出边界）
  theme = tm,        # 应用之前定义的森林图主题（控制字体、颜色、线型等外观）
  header = new_headers,  # 指定表头名称（与plot_data的列名一致，确保表头正确显示）
  col_widths = c(3, 1, 1.5, 1.5, 3, 2)  # 手动指定每列的相对宽度（6个数值对应6列，数值越大列越宽）
)

# 8. 图形美化
# 设置置信区间颜色（按方法区分）
boxcolor <- case_when(
  data$method == "Wald ratio" ~ "#E64B35",
  data$method == "Inverse variance weighted (fixed effects)" ~ "red",
  TRUE ~ "#00A087"  # 默认颜色（一般不会触发）
)

for (i in 1:nrow(data)) {
  plot <- edit_plot(plot, col = 5, row = i, which = "ci",
                    gp = gpar(fill = boxcolor[i], fontsize = 25))
}

# 加粗显著P值（<0.05）
pos_bold_pval <- which(as.numeric(gsub('<', "", data$pval)) < 0.05)
if (length(pos_bold_pval) > 0) {
  for (i in pos_bold_pval) {
    plot <- edit_plot(plot, col = 4, row = i, which = "text",
                      gp = gpar(fontface = "bold"))
  }
}

# 添加分组边框
plot <- add_border(plot, part = "header", row = 1, where = "top",
                   gp = gpar(lwd = 2))
plot <- add_border(plot, part = "header", row = lineVec,
                   gp = gpar(lwd = 1))

# 设置字体大小和居中对齐
plot_cols <- new_headers

# 调整数据行（正文）文本的字体大小
plot <- edit_plot(plot,
                  col = 1:length(plot_cols), # 针对所有列（基于新表头的列数）
                  row = 1:nrow(data),  # 针对所有数据行
                  which = "text", # 指定修改的是文本内容
                  gp = gpar(fontsize = 14)) # 将数据行文本字体大小统一设置为14pt
# 调整表头（列名）的水平对齐方式
plot <- edit_plot(plot,
                  col = 1:length(plot_cols),
                  which = "text",  # 指定修改的是表头部分（而非数据行）
                  hjust = unit(0.2, "npc"),  # 水平对齐方式（0.2表示从左侧开始20%的位置对齐，实现左偏对齐）。npc 是 "normalized parent coordinates"（标准化父坐标）的缩写，是一种相对坐标单位。
                  part = "header",
                  x = unit(0.2, "npc")) # 表头文本的x轴位置（与hjust配合确保对齐一致）

# 调整数据行（正文）文本的水平对齐方式
plot <- edit_plot(plot,
                  col = 1:length(plot_cols), # 针对所有列的数据行
                  which = "text",
                  hjust = unit(0.2, "npc"),  # 水平对齐方式（与表头保持一致的左偏对齐，视觉统一）
                  x = unit(0.2, "npc"))  # 数据行文本的x轴位置（与hjust配合确保对齐一致）

# 9. 输出森林图
output_pdf <- here(output_dir, "forest.pdf")
pdf(output_pdf, width = 25, height = 25)
print(plot)
dev.off()

message("森林图已保存至: ", output_pdf)
message("所有操作完成!")
