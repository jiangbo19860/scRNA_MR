# 加载必需包
library(GSVA)
library(clusterProfiler)
library(limma)
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)

# ------------------------------------------------------------------------------
# 0. 数据预处理（单个基因集分析也需先运行）
mg <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/Microglia_亚分4类.rds")
expr_matrix <- as.matrix(GetAssayData(mg, assay = "SCT", layer = "scale.data"))  # 表达矩阵
group_info <- mg$orig.ident  # 分组信息（Epilepsy/Healthy）
table(group_info)  # 确认分组正确

# ------------------------------------------------------------------------------
# 1. 定义当前要分析的基因集（仅修改这里即可切换其他基因集）
gmt_path <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c2.cp.reactome.v2025.1.Hs.symbols.gmt"  # 基因集路径
prefix <- "c2.reactome"  # 前缀（用于命名结果文件）

# ------------------------------------------------------------------------------
# 2. 读取并筛选基因集（针对Reactome的专属处理）
gs <- read.gmt(gmt_path)  # 读取gmt文件
gs$term <- str_remove(gs$term, "REACTOME_")  # 去除"REACTOME_"前缀，简化名称

# 转换为基因集列表，仅保留在表达矩阵中存在的基因
gs_list <- gs %>% split(.$term) %>% lapply("[[", 2)  # 按通路拆分基因
gs_list <- lapply(gs_list, function(x) intersect(x, rownames(expr_matrix)))  # 匹配表达矩阵基因
gs_list <- gs_list[sapply(gs_list, length) >= 2]  # 过滤掉基因数<2的通路
cat(paste0("有效通路数：", length(gs_list), "\n"))  # 检查是否有有效通路（至少1个）

# ------------------------------------------------------------------------------
# 3. 运行SSGSEA
ssgsea_param <- ssgseaParam(
  exprData = expr_matrix,    # 表达矩阵
  geneSets = gs_list,        # 筛选后的通路列表
  minSize = 2,               # 最小通路基因数（与筛选一致）
  maxSize = 500,             # 最大通路基因数
  alpha = 0.25,              # 尾权重指数（默认）
  normalize = TRUE,          # 标准化分数（默认）
  verbose = FALSE
)
gsva_res <- gsva(ssgsea_param)  # 运行分析
cat(paste0("SSGSEA结果维度：", dim(gsva_res), "\n"))  # 通路数×细胞数

# ------------------------------------------------------------------------------
# 4. 差异检验（Epilepsy vs Healthy）
design <- model.matrix(~0 + factor(group_info))  # 构建分组矩阵
colnames(design) <- levels(factor(group_info))  # 列名：Epilepsy、Healthy
rownames(design) <- colnames(gsva_res)  # 行名匹配细胞名

contrast <- makeContrasts(Epilepsy - Healthy, levels = design)  # 定义对比（病例-对照）
fit <- lmFit(gsva_res, design)  # 线性模型拟合
fit2 <- contrasts.fit(fit, contrast)  # 应用对比
fit2 <- eBayes(fit2)  #  empirical Bayes调整
diff_res <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")  # 提取所有差异结果

# ------------------------------------------------------------------------------
# 5. 绘制差异通路条形图
plot_df <- data.frame(
  Pathway = rownames(diff_res),
  t_value = diff_res$t,  # t值（越大表示病例组上调越显著）
  FDR = diff_res$adj.P.Val  # 校正后P值
) %>% arrange(t_value) %>% mutate(Pathway = factor(Pathway, levels = Pathway))  # 按t值排序

# 标记显著差异（FDR<0.05且|t|>2）
plot_df$signif <- ifelse(plot_df$FDR < 0.05 & abs(plot_df$t_value) > 2, "Significant", "Not")

# 绘图
p <- ggplot(plot_df, aes(x = Pathway, y = t_value, fill = signif)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +  # 横向条形图（避免通路名重叠）
  scale_fill_manual(values = c("Significant" = "#C73E1D", "Not" = "#CCCCCC")) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", alpha = 0.7) +  # 显著阈值线
  labs(
    x = "Reactome Pathways",
    y = "t-value (Epilepsy - Healthy)",
    fill = "Significant (FDR<0.05, |t|>2)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # 去除网格线
    legend.position = "bottom",    # 图例在底部
    axis.text.y = element_text(size = 8)  # 通路名字体调小（避免重叠）
  ) +
  ggtitle("Reactome Pathway Differences (Epilepsy vs Healthy Microglia)")

# ------------------------------------------------------------------------------
# 6. 保存结果
ggsave(paste0(prefix, "_SSGSEA_Barplot.png"), p, width = 12, height = 10, dpi = 300)  # 条形图
write.csv(diff_res, paste0(prefix, "_SSGSEA_Diff.csv"), row.names = TRUE)  # 差异结果表格

# 可选：将SSGSEA结果添加到Seurat对象（便于后续可视化）
mg[[paste0("GSVA_", prefix)]] <- CreateAssayObject(data = gsva_res)
saveRDS(mg, paste0("Microglia_with_", prefix, "_GSVA.rds"))

cat(paste0("\n✅ ", prefix, " 单个基因集分析完成！\n结果文件：\n",
           "- 条形图：", prefix, "_SSGSEA_Barplot.png\n",
           "- 差异通路表：", prefix, "_SSGSEA_Diff.csv\n",
           "- Seurat对象：Microglia_with_", prefix, "_GSVA.rds"))