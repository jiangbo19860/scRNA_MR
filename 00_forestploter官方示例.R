rm(list = ls())  # 清除工作空间中的所有对象
# 加载grid包（用于底层图形绘制）
library(grid)
# 加载forestploter包（用于绘制森林图）
library(forestploter)

# 从forestploter包的内置示例数据中读取CSV文件
# system.file()用于定位包内部的文件路径，extdata是包中存放示例数据的目录
dt <- read.csv(system.file("extdata", "example_data.csv", package = "forestploter"))
colnames(dt)  # 查看数据框的列名
head(dt)  # 查看数据框的前几行

# 对Subgroup（亚组）列进行缩进处理：
# 如果Placebo列的值为NA（通常表示总组，非亚组），则保持原Subgroup名称；
# 否则在亚组名称前添加3个空格，实现视觉上的层级区分（总组顶格，亚组缩进）
dt$Subgroup <- ifelse(is.na(dt$Placebo),
                      dt$Subgroup,
                      paste0("   ", dt$Subgroup))

# 将Treatment（治疗组）列中的NA值替换为空字符串，避免绘图时显示NA
dt$Treatment <- ifelse(is.na(dt$Treatment), "", dt$Treatment)
# 同理，将Placebo（安慰剂组）列中的NA值替换为空字符串
dt$Placebo <- ifelse(is.na(dt$Placebo), "", dt$Placebo)
colnames(dt)  # 查看数据框的列名
head(dt)  # 查看数据框的前几行

# 计算标准误（se）：基于log转换后的可信区间上下限和1.96（95%CI对应的Z值）
# 公式来源：se = (log(上限) - log(效应值)) / 1.96（适用于HR等对数正态分布的效应值）
dt$se <- (log(dt$hi) - log(dt$est))/1.96

# 添加一个空白列（列名为空格），用于在森林图中放置可信区间线条
# 用20个空格填充，通过调整空格数量可控制该列的宽度（影响可信区间的显示空间）
dt$` ` <- paste(rep(" ", 20), collapse = " ")

# 创建用于展示的HR（风险比）及95%CI列：
# 如果se（标准误）为NA（通常表示无效应值的行，如标题行），则显示空字符串；
# 否则用sprintf格式化字符串，保留2位小数，格式为“效应值 (下限 to 上限)”
dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$est, dt$low, dt$hi))

colnames(dt)  # 查看数据框的列名
head(dt)  # 查看数据框的前几行

# 定义森林图的主题样式
tm <- forest_theme(
  base_size = 10,  # 图形基础字体大小
  refline_gp = gpar(col = "red"),  # 参考线（ref_line）的属性：颜色设为红色
  arrow_type = "closed",  # 坐标轴两端箭头的类型：闭合箭头
  footnote_gp = gpar(col = "blue", cex = 0.6)  # 脚注的属性：蓝色，字体大小0.6倍
)

# 绘制森林图
p <- forest(
  dt[,c(1:3, 20:21)],  # 用于展示的列：第1-3列（Subgroup、Treatment、Placebo）和第20-21列（空白列、HR (95% CI)）
  est = dt$est,  # 效应值（如HR），用于绘制可信区间的中点
  lower = dt$low,  # 可信区间下限
  upper = dt$hi,  # 可信区间上限
  sizes = dt$se,  # 标准误，用于控制可信区间线条的视觉大小（非必需，影响美观）
  ci_column = 4,  # 可信区间绘制的列索引：第4列（即前面添加的空白列）
  ref_line = 1,  # 参考线的位置（HR=1为无效应值的分界线）
  arrow_lab = c("Placebo Better", "Treatment Better"),  # 坐标轴两端的箭头标签：左侧为“安慰剂更优”，右侧为“治疗组更优”
  xlim = c(0, 4),  # X轴的范围（从0到4）
  ticks_at = c(0.5, 1, 2, 3),  # X轴刻度的位置
  footnote = "This is the demo data. Please feel free to change\nanything you want.",  # 脚注文本，\n表示换行
  theme = tm  # 应用前面定义的主题样式
)

# 关闭当前所有图形设备（清除可能残留的旧绘图设备，避免干扰）
dev.off()
# 打印森林图
plot(p)
dev.off()
