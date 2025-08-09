
rm(list = ls())  # 清除工作空间
pacman::p_load(circlize, here)  # circlize包是绘制基因位置信息的环状图的包


# 设置输入文件路径
geneRT <- read.table(here("3_outputs/汇总有显著意义的MR结果/gene.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE) # 基因列表文件
geneRT <- unique(geneRT)  # 去重，确保每个基因只出现一次

genepos = read.table(here("1_data/geneREF.txt"), header = TRUE, sep = "\t", check.names = FALSE)
colnames(genepos) = c('genename', 'chr', 'start', 'end')  # 重命名列为基因名、染色体、起始位置、终止位置
genepos = genepos[, c('chr', 'start', 'end', 'genename')]  # 重新排列列的顺序
row.names(genepos) = genepos[, 'genename']  # 将基因名设为行名
genepos = genepos[as.vector(geneRT[, 1]), ]  # 根据基因列表筛选位置信息
bed0 = genepos  # 将筛选后的基因位置信息保存为bed0



# 绘制环形图并保存为PDF文件
pdf(file = here("3_outputs/circlize.pdf"), width = 6, height = 6)  # 创建保存为PDF的绘图设备

# 开始绘制环形图
circos.clear()  # 清除circlize环境
circos.initializeWithIdeogram(species = "hg38", plotType = NULL)  # 使用人类基因组版本hg38初始化环状图

# 绘制染色体背景的彩色条并注释染色体信息
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index  # 获取当前染色体的名称
  xlim = CELL_META$xlim  # 获取当前染色体的范围
  ylim = CELL_META$ylim  # 获取当前绘图的y轴范围
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(24))  # 在当前染色体区域填充随机颜色
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.6, col = "white",  # 在染色体区域中央添加染色体名称
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)  # 设置环形条高度和边框

# 绘制基因组染色体的带状图
circos.genomicIdeogram(species = "hg38", track.height = mm_h(6))  # 显示基因组的染色体图谱

# 绘制基因的位置信息并注释基因名称
circos.genomicLabels(bed0, labels.column = 4, side = "inside", cex = 0.8)  # 根据基因位置绘制标签

# 清除circlize环境
circos.clear()

# 关闭PDF输出设备
dev.off()
