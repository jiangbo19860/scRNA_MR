# 绘制人类基因组中目标基因位置-环状图（Circular Plot），2个输入：sig_genes.txt文件中的基因列表, 基因位置参考文件1_data/geneREF.txt。
# 输出1个环状图：3_outputs/circlize.pdf

rm(list = ls())  # 清除工作空间
pacman::p_load(circlize, here)  # 加载所需包

# 设置输入文件路径
geneRT <- read.table(here("3_outputs/20250810_original_diffGeneExp_p05_26genes.txt"),
                     sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
geneRT <- unique(geneRT)  # 去重，确保每个基因只出现一次

genepos <- read.table(here("1_data/geneREF.txt"),
                      header = TRUE, sep = "\t", check.names = FALSE)
colnames(genepos) <- c('genename', 'chr', 'start', 'end')  # 重命名列
genepos <- genepos[, c('chr', 'start', 'end', 'genename')]  # 调整列顺序
row.names(genepos) <- genepos[, 'genename']  # 基因名设为行名

# 根据基因列表筛选位置信息，并去除无效基因（在genepos中不存在的基因）
valid_genes <- intersect(as.vector(geneRT[, 1]), rownames(genepos))  # 保留存在位置信息的基因
genepos <- genepos[valid_genes, ]  # 筛选有效基因的位置信息
bed0 <- genepos  # 保存筛选结果

# 计算图片中实际显示的基因数（去除NA和无效条目后）
gene_count <- nrow(bed0)  # 有效基因数量
cat("图片中显示的基因数：", gene_count, "\n")

# 生成包含日期和基因数的文件名
today_date <- format(Sys.Date(), "%Y%m%d")  # 当天日期（格式：年日月）
output_filename <- paste0(today_date, "_circlize_", gene_count, "genes.pdf")  # 新文件名

# 绘制环形图并保存（使用新文件名）
pdf(file = here("3_outputs", output_filename), width = 6, height = 6)

# 开始绘制环形图
circos.clear()
circos.initializeWithIdeogram(species = "hg38", plotType = NULL)

# 绘制染色体背景及注释
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(24))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.6, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)

# 绘制染色体图谱
circos.genomicIdeogram(species = "hg38", track.height = mm_h(6))

# 绘制基因位置及名称（仅有效基因）
if (gene_count > 0) {  # 确保有基因可绘制
  circos.genomicLabels(bed0, labels.column = 4, side = "inside", cex = 0.8)
} else {
  warning("无有效基因可显示，生成空图")
}

circos.clear()
dev.off()  # 关闭PDF设备

cat("环状图已保存至：", here("3_outputs", output_filename), "\n")
