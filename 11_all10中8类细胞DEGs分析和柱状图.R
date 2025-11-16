all10 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")
head(all10)
unique(all10$cell_type_8class)
unique(all10$orig.ident)
# 1. 加载必需包（确保版本匹配，Seurat≥4.0，ggplot2≥3.4.0）
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)

# 2.3 确认分组和细胞类型的因子格式（避免分析报错）
all10$orig.ident <- factor(all10$orig.ident, levels = c("Healthy", "Epilepsy"))  # 设Healthy为对照，Epilepsy为处理
all10$cell_type_8class <- factor(all10$cell_type_8class)  # 细胞类型转为因子

# 1. 提取RNA assay的原始counts（关键：必须是"counts"槽位，非"data"）
raw_rna_counts <- GetAssayData(
  object = all10,
  assay = "RNA",    # 从RNA assay提取（SCT的输入是原始UMI）
  slot = "counts"   # "counts"=原始未标准化数据，"data"=已log归一化，不可用
)

# 2. 提取原all10的核心meta信息（保留细胞群、分组、样本ID等关键列）
# 按需保留列，至少包含cell_type_8class（细胞群）、orig.ident（分组）
meta_keep <- all10@meta.data[, c("cell_type_8class", "orig.ident", "Sample_ID")]  # 可根据你的列名调整

# 3. 重建新Seurat对象（避免原对象中7个模型的干扰）
all10_clean <- CreateSeuratObject(
  counts = raw_rna_counts,  # 原始UMI计数
  meta.data = meta_keep,    # 保留原分组信息
  project = "Epilepsy_Healthy"
)

# 验证：新对象的基因/细胞数与原all10一致（确保数据未丢失）
cat("新对象（all10_clean）：基因数", nrow(all10_clean), "，细胞数", ncol(all10_clean), "\n")
cat("原对象（all10）：基因数", nrow(all10), "，细胞数", ncol(all10), "\n")
# ✅ 正常输出：两者基因数、细胞数完全一致


# 重新执行SCT标准化（核心：生成全局统一模型）
all10_clean <- SCTransform(
  object = all10_clean,
  assay = "RNA",          # 输入：RNA assay的原始counts
  new.assay.name = "SCT", # 输出：新建SCT assay（覆盖原多模型）
  verbose = FALSE         # 关闭冗余输出，聚焦结果
)

# 关键验证：检查新对象的SCT模型数量（必须为1，否则重新执行此步）
cat("重新SCT后，模型数量：", length(all10_clean@assays$SCT@SCTModel.list), "\n")
# ✅ 正常输出：1（若仍>1，重启R后重新执行步骤1-2，避免缓存干扰）

all10_clean <- PrepSCTFindMarkers(
  object = all10_clean,
  assay = "SCT",
  verbose = FALSE
)

# 验证1：SCT模型数量（必须=1）
if (length(all10_clean@assays$SCT@SCTModel.list) == 1) {
  cat("✅ 多模型问题彻底解决！当前模型数量：1\n")
} else {
  stop("❌ 模型数量仍>1，请重启R后重新执行步骤1-2")
}

# 验证2：细胞群和分组信息是否保留（确保后续分析可用）
cat("\n细胞群类型：", paste(unique(all10_clean$cell_type_8class), collapse = ", "), "\n")
cat("分组类型：", paste(unique(all10_clean$orig.ident), collapse = ", "), "\n")
# ✅ 正常输出：8个细胞群、2个分组（与原all10一致）

# saveRDS(all10_clean, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/1116_all10重新SCT标准化后.rds")


all10_clean <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/1116_all10重新SCT标准化后.rds")

# 重新定义diff_params
diff_params <- list(
  test.use = "wilcox",        # Wilcoxon秩和检验
  only.pos = FALSE,          # 保留上下调基因
  min.pct = 0.25,            # 至少25%细胞表达
  logfc.threshold = 0.5,     # log2FC绝对值>0.5
  assay = "SCT",             # 基于SCT assay
  group.by = "orig.ident",   # 分组列
  ident.1 = "Epilepsy",      # 处理组
  ident.2 = "Healthy",       # 对照组
  verbose = FALSE            # 关闭冗余输出（新增，避免刷屏）
)

# 重新执行循环分析（核心：手动指定每个参数）
cell_types <- unique(all10_clean$cell_type_8class)
deg_list <- list()

for (ct in cell_types) {
  # 1. 子集化当前细胞群
  ct_subset <- subset(all10_clean, subset = cell_type_8class == ct)
  
  # 2. 过滤分组细胞数不足的情况（每组至少5个细胞）
  healthy_num <- sum(ct_subset$orig.ident == "Healthy")
  epilepsy_num <- sum(ct_subset$orig.ident == "Epilepsy")
  if (healthy_num < 5 | epilepsy_num < 5) {
    message(paste0("⚠️  跳过细胞群：", ct, "（Healthy组", healthy_num, "个，Epilepsy组", epilepsy_num, "个，需各≥5个）"))
    next
  }
  
  # 3. 执行差异分析（关键：手动传入每个参数，无语法冲突）
  deg <- FindMarkers(
    object = ct_subset,
    test.use = diff_params$test.use,        # 逐个指定参数
    only.pos = diff_params$only.pos,        # 对应diff_params的键
    min.pct = diff_params$min.pct,          # 确保参数名完全匹配
    logfc.threshold = diff_params$logfc.threshold,
    assay = diff_params$assay,
    group.by = diff_params$group.by,
    ident.1 = diff_params$ident.1,
    ident.2 = diff_params$ident.2,
    verbose = diff_params$verbose
  )
  
  # 4. 筛选显著DEGs（用基础R替代rownames_to_column，无dplyr依赖）
  deg_signif <- data.frame(
    gene = rownames(deg),  # 手动将行名（基因名）转为第一列“gene”
    deg,                   # 原差异分析结果（含avg_log2FC、p_val_adj等）
    check.names = FALSE    # 避免列名自动修改（如去除特殊字符）
  )
  
  # 后续筛选逻辑不变（过滤p_val_adj < 0.05 + 定义差异方向）
  deg_signif <- deg_signif %>%
    filter(p_val_adj < 0.05) %>%  # 校正后P值<0.05
    mutate(
      cell_type = ct,
      DEG_direction = case_when(
        avg_log2FC > 0.5 ~ "Up",    # 上调：Epilepsy > Healthy
        avg_log2FC < -0.5 ~ "Down", # 下调：Epilepsy < Healthy
        TRUE ~ "NS"                 # 不显著（理论上已过滤）
      )
    ) %>%
    filter(DEG_direction != "NS")   # 移除不显著基因
  
  # （可选）清除行名（避免后续混淆）
  rownames(deg_signif) <- NULL
  
  # 5. 存储结果并打印进度
  deg_list[[ct]] <- deg_signif
  message(paste0("✅ 完成细胞群：", ct, "（上调基因", sum(deg_signif$DEG_direction=="Up"), "个，下调基因", sum(deg_signif$DEG_direction=="Down"), "个）"))
}

# 3. 可视化：水平条形图（Horizontal Bar Chart），
deg_all <- bind_rows(deg_list)
deg_count <- deg_all %>%
  group_by(cell_type, DEG_direction) %>%
  summarise(gene_num = n(), .groups = "drop") %>%
  complete(cell_type, DEG_direction, fill = list(gene_num = 0))

# 调整deg_count数据：Down的gene_num设为负，Up为正，便于水平条形图左右区分
# 先恢复deg_count的原始数值（将Down的负数转回正数，方便计算总数）
deg_count <- deg_count %>%
  mutate(gene_num = abs(gene_num))

# 计算每个细胞类型的总DEGs数（Up + Down）
cell_type_total <- deg_count %>%
  group_by(cell_type) %>%
  summarise(total_degs = sum(gene_num), .groups = "drop")

# 按总DEGs数降序排序细胞类型
sorted_cell_types <- cell_type_total %>%
  arrange(desc(total_degs)) %>%
  pull(cell_type)

# 将cell_type转换为因子，levels设为排序后的顺序
deg_count <- deg_count %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(sorted_cell_types)),
    gene_num = ifelse(DEG_direction == "Down", -gene_num, gene_num)  # 再将Down设为负数，用于水平条形图
  )

# 计算差异基因数的最大绝对值并向上取整到100的倍数
max_down <- max(deg_count[deg_count$DEG_direction == "Down", "gene_num"] * -1)
max_up <- max(deg_count[deg_count$DEG_direction == "Up", "gene_num"])
max_abs <- max(max_down, max_up)
max_rounded <- ceiling(max_abs / 100) * 100  # 向上取整到最近的100倍数

p <- ggplot(deg_count, aes(x = gene_num, y = cell_type, fill = DEG_direction)) +
  geom_col(width = 0.7, color = "black", size = 0.1, lineend = "square") +
  scale_fill_manual(values = c("Down" = "#0072B2", "Up" = "#D55E00")) +
  scale_x_continuous(
    breaks = seq(-max_rounded, max_rounded, by = 100),  # 每100一个刻度
    labels = abs(seq(-max_rounded, max_rounded, by = 100)),  # 刻度标签显示绝对值
    limits = c(-max_rounded, max_rounded)  # x轴范围按100倍数取整
  ) +
  labs(x = "DEGs", y = "Cell Type", title = "Epilepsy vs Healthy") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "plain", size = 10, color = "black"),
    axis.text.y = element_text(face = "plain", size = 10, color = "black"),
    axis.title.x = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
    plot.margin = margin(20, 20, 20, 20, unit = "pt")
  )

p

ggsave("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_8群细胞DEGs.pdf", p, width = 10, height = 6, bg = "white")

# 1. 保存所有细胞群合并后的差异基因总表（路径去除多余空格，避免保存失败）
write.csv(
  x = deg_all,
  file = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_所有细胞群DEGs总表.csv",
  row.names = FALSE,
  quote = FALSE
)

# 2. 按细胞群单独保存差异基因表
for (ct in names(deg_list)) {
  write.csv(
    x = deg_list[[ct]],
    file = paste0("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/1116_", ct, "_DEGs.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

# 保存完成提示
message("✅ 差异基因CSV文件保存完成！")
message("📁 合并总表：1116_所有细胞群DEGs总表.csv")
message("📁 分群表：每个细胞群对应独立CSV（如1116_Excitatory_Neurons_DEGs.csv）")

