# 关联单细胞数据的标志基因与 TCGA 肿瘤数据，分析不同聚类的标志基因在 LUAD（肺腺癌）和 LUSC（肺鳞癌）中的表达差异。整体流程为：整理单细胞数据的细胞类型注释（合并相似聚类）→ 2. 提取每个聚类的核心标志基因 → 3. 分析这些标志基因在 LUAD 和 LUSC 肿瘤中的平均表达差异（计算倍数变化和 p 值）→ 4. 绘制散点图展示差异，通过对角线判断表达高低（LUSC vs LUAD）。

# TCGA 的基因表达数据来自肿瘤组织块的整体测序（bulk RNA-seq），即从一块肿瘤组织中提取所有细胞的总 RNA 进行测序，得到的是该组织中所有细胞（包括癌细胞、免疫细胞、基质细胞等）的基因表达平均值。TCGA 的 bulk 数据反映肿瘤整体的分子特征（如某基因在某类癌症中普遍高表达），且样本量大、临床信息完整，适合关联预后、药物反应等宏观表型。单细胞数据可解析bulk 数据中 “平均信号” 的细胞来源（如某基因的高表达是来自癌细胞还是免疫细胞），揭示微观层面的机制。

rm(list = ls())

# 加载所需R包：Seurat用于单细胞数据分析；tidyverse用于数据处理和可视化
library(Seurat)
library(tidyverse)

# 加载单细胞RNA-seq数据对象scRNA（包含基因表达和细胞元信息）
load("scRNA.Rds")

# 从scRNA中提取元信息（metadata），每行代表一个细胞，包含聚类、样本等信息
metadata <- scRNA@meta.data

# 在metadata中新增"cell"列，初始值为聚类名称（ClusterName，如"T cells 1"等）
metadata$cell <- metadata$ClusterName

# 对"cell"列进行细胞类型合并（将相似聚类合并为更广泛的细胞类型）：
# 1. 包含"T cells"或"natural killer"的聚类→统一标记为"T cells"（T细胞及NK细胞）
metadata[grep("T cells|natural killer", metadata$ClusterName), ]$cell <- "T cells"
# 2. 包含"B cells"或"granulocytes"的聚类→统一标记为"B cells"（B细胞及粒细胞）
metadata[grep("B cells|granulocytes", metadata$ClusterName), ]$cell <- "B cells"
# 3. 包含"fibroblasts"的聚类→统一标记为"Fibroblasts"（成纤维细胞）
metadata[grep("fibroblasts", metadata$ClusterName), ]$cell <- "Fibroblasts"
# 4. 包含"endothelial cell"的聚类→统一标记为"Endothelial"（内皮细胞）
metadata[grep("endothelial cell", metadata$ClusterName), ]$cell <- "Endothelial"
# 5. 包含"epithelial cell|EC|basal cells"的聚类→统一标记为"Epithelial"（上皮细胞）
metadata[grep("epithelial cell|EC|basal cells", metadata$ClusterName), ]$cell <- "Epithelial"
# 6. 包含"cancer cells"的聚类→统一标记为"Cancer"（癌细胞）
metadata[grep("cancer cells", metadata$ClusterName), ]$cell <- "Cancer"
# 7. 包含"alveolar"的聚类→统一标记为"Alveolar"（肺泡细胞）
metadata[grep("alveolar", metadata$ClusterName), ]$cell <- "Alveolar"
# 8. 包含"macrophages|Langerhans|mast cells|dendritic"的聚类→统一标记为"Myleoid"（髓系细胞，如巨噬细胞、树突状细胞等）
metadata[grep("macrophages|Langerhans|mast cells|dendritic", metadata$ClusterName), ]$cell <- "Myleoid"
# 9. 包含"erythroblasts|secretory club cells"的聚类→统一标记为"basel"（基底细胞相关）
metadata[grep("erythroblasts|secretory club cells", metadata$ClusterName), ]$cell <- "basel"

# 将整理后的metadata放回scRNA对象中（更新细胞类型注释）
metadata -> scRNA@meta.data

# 设置scRNA的默认身份（Idents）为"cell"列（即合并后的细胞类型）
Idents(scRNA) <- "cell"

# 绘制降维图（默认UMAP或t-SNE），展示不同细胞类型的分布（按"cell"列着色）
DimPlot(scRNA)

# 注释：寻找所有细胞类型的标志基因（仅保留上调基因，log2FC>0.5，最小表达比例0.25）
# 由于可能已提前运行并保存结果，这里注释掉
# sc.marker <- FindAllMarkers(scRNA,
#                            only.pos = T,
#                            logfc.threshold = 0.5,
#                            min.pct = 0.25)

# 加载预保存的成纤维细胞相关单细胞数据（可能是scRNA的子集，聚焦于成纤维细胞）
load("scRNA_fibro.Rds")

# 重新计算标志基因：寻找每个聚类的标志基因（仅上调基因，log2FC>0.5，最小表达比例0.25）
sc.marker <- FindAllMarkers(scRNA, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)

# 对标志基因进行筛选：按聚类分组，每个聚类仅保留log2FC最高的5个基因（精简结果）
sc.marker <- sc.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# 加载TCGA数据（包含LUAD、LUSC等肿瘤的表达数据，如LUAD_tumor、LUSC_tumor等）
load("TCGA.Rds")

# 创建空数据框，用于存储每个聚类的分析结果（LUAD和LUSC的平均表达、倍数变化、p值等）
data_plot_all <- data.frame(
  LUAD = NULL,       # LUAD肿瘤中标志基因的平均表达
  LUSC = NULL,       # LUSC肿瘤中标志基因的平均表达
  fc = NULL,         # 倍数变化（LUSC平均 / LUAD平均）
  p = NULL,          # 差异显著性p值（Wilcoxon检验）
  cluster = NULL     # 聚类编号
)

# 循环分析每个聚类（0到5号聚类，共6个聚类）
for (x in 0:5) {
  # 提取当前聚类（x）的标志基因
  sc.marker.tmp <- subset(sc.marker, cluster == x)
  
  # 筛选出同时存在于lung（可能是正常肺组织表达矩阵）和当前聚类标志基因中的基因
  # 目的：确保标志基因在TCGA数据中存在表达记录
  marker <- intersect(rownames(lung), sc.marker.tmp$gene)
  
  # 提取LUAD肿瘤中这些标志基因的表达数据，并转换为数据框
  data_LUAD <- as.data.frame(LUAD_tumor[marker, ])
  
  # 计算每个LUAD样本中，当前聚类标志基因的平均表达量（新增"avr"行存储平均值）
  for (i in 1:ncol(data_LUAD)) {  # i为样本索引
    data_LUAD["avr", i] <- mean(as.numeric(data_LUAD[1:length(marker), i]))  # 对该样本的所有标志基因取平均
  }
  
  # 提取LUSC肿瘤中这些标志基因的表达数据，并转换为数据框
  data_LUSC <- as.data.frame(LUSC_tumor[marker, ])
  
  # 计算每个LUSC样本中，当前聚类标志基因的平均表达量（新增"avr"行存储平均值）
  for (i in 1:ncol(data_LUSC)) {
    data_LUSC["avr", i] <- mean(as.numeric(data_LUSC[1:length(marker), i]))
  }
  
  # 提取LUSC样本的平均表达值（"avr"行的所有样本）
  a <- as.numeric(data_LUSC["avr", ])
  # 提取LUAD样本的平均表达值
  b <- as.numeric(data_LUAD["avr", ])
  
  # 用Wilcoxon秩和检验比较LUSC和LUAD中标志基因平均表达的差异，获取p值
  p <- wilcox.test(a, b, na.rm = T)$p.value
  # 计算倍数变化（LUSC平均表达 / LUAD平均表达）
  fc <- mean(a) / mean(b)
  
  # 将当前聚类的结果整理为数据框
  data_plot <- data.frame(
    LUSC = mean(a),    # LUSC中标志基因的平均表达（所有样本的均值）
    LUAD = mean(b),    # LUAD中标志基因的平均表达（所有样本的均值）
    fc = fc,           # 倍数变化
    p = p,             # 差异p值
    cluster = x        # 聚类编号
  )
  
  # 将当前聚类结果合并到总结果数据框中
  data_plot_all <- rbind(data_plot, data_plot_all)
}

# 对p值进行多重检验校正（FDR方法），新增padj列（校正后的p值）
data_plot_all$padj <- p.adjust(data_plot_all$p, method = "fdr")

# 设置散点图的坐标轴范围：比LUAD和LUSC的最小/最大值略宽（上下各加0.05），确保所有点都能显示
scale_lim <- c(
  (min(data_plot_all[, c("LUAD", "LUSC")]) - 0.05),  # 最小值
  (max(data_plot_all[, c("LUAD", "LUSC")]) + 0.05)   # 最大值
)

# 加载ggrepel包，用于避免散点图中标签重叠
library(ggrepel)

# 对倍数变化进行分组标记：fc>1→">1"（LUSC中表达更高）；fc<1→"<1"（LUAD中表达更高）
data_plot_all$FC <- ifelse(data_plot_all$fc > 1, ">1", "<1")

# 为每个聚类添加标签（如"C0"、"C1"等）
data_plot_all$label <- paste0("C", data_plot_all$cluster)

# 绘制散点图，展示每个聚类的标志基因在LUAD和LUSC中的表达差异
ggplot(data = data_plot_all) +
  # 绘制散点：x轴=LUSC平均表达，y轴=LUAD平均表达，颜色按FC分组（">1"为红色，"<1"为黑色）
  geom_point(aes(x = LUSC, y = LUAD, color = FC), alpha = 0.5) +
  # 自定义颜色：黑色（fc<1）和红色（fc>1）
  scale_color_manual(values = c("black", "red")) +
  # 添加标签：用聚类标签（如"C0"）标注每个点，避免重叠
  geom_text_repel(aes(x = LUSC, y = LUAD, label = label)) +
  # 使用白色背景主题
  theme_bw() +
  # 设置x轴和y轴范围（统一范围，便于比较）
  xlim(scale_lim) + ylim(scale_lim) +
  # 添加对角线（斜率=1，截距=0）：线上的点表示LUAD和LUSC表达相等；线以上的点表示LUAD表达更高，线以下表示LUSC表达更高
  geom_abline(slope = 1, intercept = 0, lty = 2)  # 虚线

# 该散点图的核心价值是 “桥接单细胞亚群与临床肿瘤亚型”：通过标志基因的表达差异，揭示哪些细胞亚群与 LUAD 或 LUSC 的分子特征更相关，为后续研究（如亚型特异性治疗靶点开发、肿瘤微环境差异机制）提供方向。例如，若 LUSC 高表达的 “C2” 聚类是促癌成纤维细胞亚群，则可进一步研究该亚群在 LUSC 中的功能及干预策略。
# 1. 坐标轴与对角线的含义
# X 轴：表示该聚类标志基因在 LUSC 肿瘤 中的平均表达量（所有 LUSC 样本的均值）。
# Y 轴：表示该聚类标志基因在 LUAD 肿瘤 中的平均表达量（所有 LUAD 样本的均值）。
# 对角线（虚线，斜率 = 1）：
# 线上的点：表示该聚类标志基因在 LUAD 和 LUSC 中的表达量 相等（LUSC 平均 = LUAD 平均）。
# 线 上方 的点：表示标志基因在 LUAD 中表达更高（LUAD 平均 > LUSC 平均）。
# 线 下方 的点：表示标志基因在 LUSC 中表达更高（LUSC 平均 > LUAD 平均）。
# 2. 点的颜色与标签
# 颜色：
# 红色点：FC > 1（LUSC 平均表达 / LUAD 平均表达 > 1），即标志基因在 LUSC 中表达更高。
# 黑色点：FC < 1，即标志基因在 LUAD 中表达更高。
# 标签（如 “C0”“C1”）：对应单细胞数据中的聚类编号（0~5），代表不同的细胞亚群（如某类 T 细胞亚群、成纤维细胞亚群等）。
# 3. 核心解读逻辑
# 通过点的位置和颜色，可直接判断 “哪个细胞亚群的标志基因在两种肺癌亚型中存在显著表达差异”，具体如下：
# （1）红色点（LUSC 高表达）
# 若某聚类（如 “C2”）为红色点且位于对角线下方，说明该亚群的标志基因在 LUSC 中表达显著高于 LUAD。
# 生物学意义：该细胞亚群可能在 LUSC 的肿瘤微环境中更活跃或更富集（如 LUSC 中某类促癌成纤维细胞亚群更活跃）。
# （2）黑色点（LUAD 高表达）
# 若某聚类（如 “C0”）为黑色点且位于对角线上方，说明该亚群的标志基因在 LUAD 中表达显著高于 LUSC。
# 生物学意义：该细胞亚群可能与 LUAD 的特性更相关（如 LUAD 中某类免疫抑制性 T 细胞亚群更富集）。
# （3）靠近对角线的点
# 若点接近对角线（无论颜色），说明该亚群的标志基因在两种亚型中表达差异较小，可能是两种肺癌共有的细胞亚群（如通用的巨噬细胞亚群）。
# 4. 结合临床与单细胞数据的解读
# 若某聚类（如 “C3”）为红色点且远离对角线，提示该亚群的标志基因是 LUSC 的 “特征性分子标签”，可能参与 LUSC 的发生发展（如驱动 LUSC 的侵袭转移）。
# 若某聚类（如 “C1”）为黑色点且差异显著，可能与 LUAD 的特定表型（如靶向药耐药、免疫逃逸）相关，可进一步探索其作为 LUAD biomarkers 的潜力。
