rm(list = ls())
pacman::p_load(
  here, Seurat, SCopeLoomR, loomR, dplyr, ggplot2, cowplot, patchwork, ggrepel, SingleR, celldex, BiocParallel, scater, scran, SingleCellExperiment, data.table, hdf5r, Matrix, tidyr, stringr, forcats, pheatmap, RColorBrewer, viridis, ggridges, tibble
)
here()

# 0. 定义路径，加载数据 ------------------------------------------------------------
rds_file <- here("1_data", "E-MTAB-6149", "scRNA_processed.Rds")

# 定义输出路径函数，所有输出文件保存在3_outputs目录下的E-MTAB-6149_outputs子目录中。
output_path <- function(...) {
  here("3_outputs", "E-MTAB-6149_outputs", ...)
}

# 创建输出目录（如果不存在）
if (!dir.exists(output_path())) {
  dir.create(output_path(), recursive = TRUE)
  message("已创建输出目录:", output_path())
}

# 读取Rds文件
load(rds_file) # 读取Rds文件，加载数据到当前工作空间
ls()

metadata <- scRNA@meta.data

# Fig1b_t-SNE降维图，按CellFromTumor着色 -----------------------------------------------------------
unique(scRNA$CellFromTumor)
DimPlot(scRNA, reduction = "tsne", group.by = "CellFromTumor") +
  # 根据unique(scRNA$CellFromTumor)返回的实际值 "1" 和 "0" 调整颜色映射
  scale_color_manual(
    values = c("1" = "#2D5474", "0" = "#72C667"), # "1" 对应肿瘤，"0" 对应非肿瘤
    labels = c("Tumor", "Non-malignant"), # 自定义图例标签（顺序与 values 一致）
    limits = c("1", "0") # 控制图例顺序（Tumor 在上）
  ) +
  # 将图例形状改为方块（通过guides设置）
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) + # 15=方块，16=圆，17=三角.
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.position = c(.01, .05), # 图例位置（左下角），(0,0)为左下角，(1,1)为右上角
    legend.key = element_rect(fill = NA) # 去除图例背景
  ) +
  labs(title = "Sample Origin")

ggsave(filename = output_path("Fig1b_1_tSNE_CellFromTumor.pdf"), height = 8, width = 8)

# Fig1b中的第2个图：按患者编号绘制t-SNE图 --------------------------------------------------------------------
unique(scRNA@meta.data$PatientNumber) # "1" "2" "3" "4" "5"
DimPlot(scRNA, reduction = "tsne", group.by = "PatientNumber") +
  # 根据unique(scRNA$CellFromTumor)返回的实际值 "1" 和 "0" 调整颜色映射
  scale_color_manual(
    values = c("1" = "#268A24", "2" = "#8EEE8B", "3" = "#FB6346", "4" = "#FBD51A", "5" = "#28507D"),
    labels = c("1", "2", "3", "4", "5"), # 自定义图例标签（顺序与values一致）
    # limits = c("1", "0") # 控制图例顺序（Tumor 在上）
  ) +
  # 将图例形状改为方块（通过guides设置）
  guides(color = guide_legend(
    override.aes = list(shape = 15, size = 5), # 15=方块，16=圆，17=三角
    nrow = 2,
    byrow = TRUE,
    # title = "Patient"
  )) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.position = c(.01, .05), # 图例位置（左下角），(0,0)为左下角，(1,1)为右上角
    legend.key = element_rect(fill = NA) # 去除图例背景
  ) +
  labs(title = "Patient")

ggsave(filename = output_path("Fig1b_2_tSNE_PatientNumber.pdf"), height = 8, width = 8)

# Fig1b中的第3个图：按细胞类型dbCluster绘制t-SNE图 --------------------------------------------------------------------
unique(metadata$dbCluster) # "B cell" "CD4 T cell" "CD8 T cell" "Dendritic cell" "Endothelial cell" "Fibroblast" "Macrophage")

metadata$cell <- metadata$ClusterName

metadata[grep("T cells|natural killer", metadata$ClusterName), ]$cell <- "T cells"
metadata[grep("B cells|granulocytes", metadata$ClusterName), ]$cell <- "B cells"
metadata[grep("fibroblasts", metadata$ClusterName), ]$cell <- "Fibroblasts"
metadata[grep("endothelial cell", metadata$ClusterName), ]$cell <- "Endothelial"
metadata[grep("epithelial cell|EC|basal cells", metadata$ClusterName), ]$cell <- "Epithelial"
metadata[grep("cancer cells", metadata$ClusterName), ]$cell <- "Cancer"
metadata[grep("alveolar", metadata$ClusterName), ]$cell <- "Alveolar"
metadata[grep("macrophages|Langerhans|mast cells|dendritic", metadata$ClusterName), ]$cell <- "Myleoid"
metadata[grep("erythroblasts|secretory club cells", metadata$ClusterName), ]$cell <- NA

table(metadata$cell) # 每个cell类型对应的细胞数
scRNA$cell <- metadata$cell

DimPlot(scRNA,
  reduction = "tsne", # 降维方法为t-SNE
  group.by = "cell", label = TRUE # 按新注释的"cell"列（合并后的细胞类型）分组着色
) +
  scale_color_manual(values = c(
    "#3A6135",
    "#E7A649",
    "#4C6679",
    "#814ECC",
    "#6466A5",
    "#A14299",
    "#F58C70",
    "#B4434E"
  )) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
    legend.position = "none"
  ) +
  labs(title = "Cell type")

ggsave(filename = output_path("Fig1b3tSNE_细胞类型.pdf"), height = 8, width = 8)

# Fig1b中的第4个图：按转录本计数nUMI绘制t-SNE图 --------------------------------------------------------------------
FeaturePlot(scRNA,
  reduction = "tsne",
  features = "nUMI"
) +
  scale_color_gradient(high = "darkblue", low = "gray") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    legend.position = "none"
  ) +
  labs(title = "Transcript count")

ggsave(filename = output_path("Fig1b4_UMIs.pdf"), height = 8, width = 8)

# 在一个 PDF 文件中批量生成并合并 8 个基因的表达量 t-SNE 图（每个基因一张图，按 4 列 2 行排列），用于展示特定基因在单细胞群体中的表达分布。

# pdf() 会创建一个 “临时输出流”，所有后续绘图操作都会写入这个流，直到用 dev.off() 关闭，否则 PDF 文件可能损坏或无法正常打开。
pdf(output_path("Fig1c.pdf"), height = 6, width = 12)

p <- FeaturePlot(scRNA,
  reduction = "tsne", # 指定降维方法为t-SNE（使用之前计算的t-SNE坐标）
  features = c( # 要展示的基因列表（共8个基因）
    "CLDN18",
    "CLDN5",
    "CAPS",
    "COL1A1",
    "CD79A",
    "LYZ",
    "CD3D",
    "EPCAM"
  ),
  ncol = 4, # 图的排列方式：按4列排列（8个基因会自动排成2行4列）
  cols = c("gray", "red") # 颜色映射：灰色表示低表达，红色表示高表达
)
p & theme(
  panel.border = element_rect( # 为每张图添加黑色边框
    fill = NA, # 边框内不填充颜色（透明）
    color = "black", # 边框颜色为黑色
    linewidth = 1, # 边框线宽
    linetype = "solid"
  ), # 边框为实线
  legend.position = "none" # 不显示图例
)
dev.off() # 关闭 PDF 设备，保存图形到文件
