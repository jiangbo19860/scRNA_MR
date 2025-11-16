rm(list = ls())
gc()

# 饼图和柱状图可视化细胞占比 -----------------------------------------------------------
# 完整代码：细胞类型比例可视化（柱状图+饼图，按orig.ident和Sample_ID分组）
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)  # 用于颜色渐变
library(forcats)

all10 <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")

DefaultAssay(all10)

# head(all10@meta.data)
# unique(all10$cell_type_8class)
# all10$Barcode     <- NULL
# all10$Sample_Name <- NULL
# all10$Donor.ID    <- NULL
# head(all10)
# saveRDS(all10, "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113/00_processed_data/5_combined_10_去掉簇9和簇12后.rds")
# --------------------------
# 第一步：数据预处理（提取细胞类型、orig.ident、Sample_ID的对应关系）
# --------------------------
# 从Seurat元数据中提取核心信息（排除Cell.type为NA的细胞）
meta_plot <- all10@meta.data %>%
  select(orig.ident, Sample_ID, cell_type_8class) %>%  
  filter(!is.na(cell_type_8class)) %>%                 # 排除未注释细胞
  mutate(
    # 统一细胞类型名称（确保与之前一致，可选）
    cell_type_8class = factor(
      cell_type_8class,
      levels = c("Excitatory_Neurons", "Inhibitory_Neurons", "Oligodendrocytes", 
                 "Astrocytes", "Microglia", "OPCs", "Endothelial_Cells", "Pericytes")
    )
  )

# 定义配色方案（8种细胞类型，颜色区分度高）
cell_type_8class_colors <- c(
  "Excitatory_Neurons" = "#E74C3C",    # 红色：兴奋性神经元
  "Inhibitory_Neurons" = "#F39C12",    # 橙色：抑制性神经元
  "Oligodendrocytes" = "#27AE60",       # 绿色：少突胶质细胞
  "Astrocytes" = "#3498DB",            # 蓝色：星形胶质细胞
  "Microglia" = "#9B59B6",             # 紫色：小胶质细胞
  "OPCs" = "#1ABC9C",                  # 青色：OPCs
  "Endothelial_Cells" = "#F1C40F",           # 黄色：内皮细胞
  "Pericytes" = "#95A5A6"              # 灰色：周细胞
)

# --------------------------
# 第二步：按orig.ident（健康/疾病）分组的可视化
# --------------------------
# 1. 计算orig.ident分组下的细胞类型数量及占比
orig_count <- meta_plot %>%
  group_by(orig.ident,cell_type_8class) %>%
  summarise(cell_num = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(
    cell_pct = cell_num / sum(cell_num) * 100,  # 计算各类型占比（%）
    pct_label = sprintf("%.1f%%\n(n=%d)", cell_pct, cell_num)  # 占比+数量标签
  )

# 2. 按orig.ident分组的柱状图（堆叠+百分比）
p_orig_bar <- orig_count %>%
  ggplot(aes(x = orig.ident, y = cell_pct, fill =cell_type_8class)) +
  geom_col(position = "stack", alpha = 0.9) +
  geom_text(
    aes(label = ifelse(cell_pct >= 5, pct_label, "")),  # 仅显示占比≥5%的标签（避免拥挤）
    position = position_stack(vjust = 0.5),
    size = 3, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = cell_type_8class_colors) +
  labs(
    title = "Cell Type Proportion by Group (orig.ident)",
    x = "Group (Healthy/Epilepsy)",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )
p_orig_bar
dev.off()

# 3. 按orig.ident分组的饼图（每个组一个饼图）
p_orig_pie <- orig_count %>%
  ggplot(aes(x = "", y = cell_pct, fill =cell_type_8class)) +
  geom_col(position = "stack", width = 1, alpha = 0.9) +
  geom_text(
    aes(label = ifelse(cell_pct >= 5, pct_label, "")),
    position = position_stack(vjust = 0.5),
    size = 3, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = cell_type_8class_colors) +
  facet_wrap(~orig.ident, ncol = 2) +  # 按组拆分饼图
  coord_polar("y") +  # 转为饼图
  labs(
    title = "Cell Type Composition by Group (orig.ident)",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_blank(),        # 隐藏饼图的轴文本
    axis.title = element_blank(),       # 隐藏轴标题
    panel.grid = element_blank(),       # 隐藏网格线
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")  # 组名称加粗
  )

p_orig_pie
dev.off()
# --------------------------
# 第三步：按Sample_ID（每个样本）分组的可视化
# --------------------------
# 1. 计算Sample_ID分组下的细胞类型数量及占比
library(dplyr)
library(ggplot2)
library(forcats)  # 用于因子操作

sample_count <- meta_plot %>%
  group_by(Sample_ID, cell_type_8class) %>%
  summarise(cell_num = n(), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(
    cell_pct  = cell_num / sum(cell_num) * 100,
    pct_label = sprintf("%.1f%%", cell_pct)
  ) %>%
  ungroup() %>%
  mutate(
    # 先把 Sample_ID 变成因子，再反转顺序
    Sample_ID = fct_rev(factor(Sample_ID))
  )

p_sample_bar <- sample_count %>%
  ggplot(aes(x = Sample_ID, y = cell_pct, fill = cell_type_8class)) +
  geom_col(position = "stack", alpha = 0.9) +
  geom_text(
    aes(label = ifelse(cell_pct >= 8, pct_label, "")),
    position = position_stack(vjust = 0.5),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = cell_type_8class_colors) +
  labs(
    title = "Cell Type Proportion by Sample (Sample_ID)",
    x = "Sample ID",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 9),
    legend.position = "right",
    legend.text  = element_text(size = 10)
  )

p_sample_bar
dev.off()


# 3. 按Sample_ID分组的饼图（每个样本一个饼图，小尺寸排列）
p_sample_pie <- sample_count %>%
  ggplot(aes(x = "", y = cell_pct, fill =cell_type_8class)) +
  geom_col(position = "stack", width = 1, alpha = 0.9) +
  geom_text(
    aes(label = ifelse(cell_pct >= 8, pct_label, "")),
    position = position_stack(vjust = 0.5),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = cell_type_8class_colors) +
  facet_wrap(~Sample_ID, ncol = 4) +  # 4列排列，根据样本数量调整ncol
  coord_polar("y") +
  labs(
    title = "Cell Type Composition by Sample (Sample_ID)",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 8, face = "bold")  # 样本名小字，避免重叠
  )
p_sample_pie

dev.off()
# --------------------------
# 第四步：保存所有图表（高分辨率，便于论文/报告使用）
# --------------------------
# 定义保存路径（与之前分析路径一致）
save_path <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/8类细胞类型注释后完整可视化/"

# 1. 保存orig.ident柱状图
ggsave(
  filename = paste0(save_path, "cell_type_8class_proportion_orig_bar.png"),
  plot = p_orig_bar,
  width = 10, height = 6, dpi = 300
)

# 2. 保存orig.ident饼图
ggsave(
  filename = paste0(save_path, "cell_type_8class_composition_orig_pie.png"),
  plot = p_orig_pie,
  width = 12, height = 6, dpi = 300
)

# 3. 保存Sample_ID柱状图
ggsave(
  filename = paste0(save_path, "cell_type_8class_proportion_sample_bar.png"),
  plot = p_sample_bar,
  width = 12, height = 8, dpi = 300  # 横向图需更高高度
)

# 4. 保存Sample_ID饼图
ggsave(
  filename = paste0(save_path, "cell_type_8class_composition_sample_pie.png"),
  plot = p_sample_pie,
  width = 16, height = 12, dpi = 300  # 多样本饼图需更大尺寸
)

# --------------------------
# 第五步：显示图表（可选，在RStudio中预览）
# --------------------------
print(p_orig_bar)
print(p_orig_pie)
print(p_sample_bar)
print(p_sample_pie)

cat("\n✅ 所有细胞类型比例图已保存至：", save_path, "\n")
cat("包含4类图表：\n")
cat("1. 按orig.ident分组的堆叠柱状图\n")
cat("2. 按orig.ident分组的饼图\n")
cat("3. 按Sample_ID分组的横向柱状图\n")
cat("4. 按Sample_ID分组的多饼图\n")

