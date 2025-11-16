# 读取→筛选→SSGSEA→差异检验→可视化→保存

# 加载必需包
library(GSVA)
library(clusterProfiler)
library(limma)
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)

# ------------------------------------------------------------------------------
# 0. 数据预处理（仅运行一次）
mg <- readRDS("/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1113重新分析/1113outputs/1115/Microglia_亚分4类.rds")
expr_matrix <- as.matrix(GetAssayData(mg, assay = "SCT", layer = "scale.data"))  # 表达矩阵
group_info <- mg$orig.ident  # 分组信息（Epilepsy/Healthy）

# 定义所有基因集文件路径（按顺序排列）
gmt_files <- c(
  "c2.all" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c2.all.v2025.1.Hs.symbols.gmt",
  "c2.kegg_legacy" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
  "c2.kegg_medicus" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt",
  "c2.reactome" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
  "c5.all" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c5.all.v2025.1.Hs.symbols.gmt",
  "c7.all" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/c7.all.v2025.1.Hs.symbols.gmt",
  "h.all" = "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/3_Refdatas/h.all.v2025.1.Hs.symbols.gmt"
)

# ------------------------------------------------------------------------------
# 1. 封装SSGSEA批量分析函数（核心简化步骤）
run_ssgsea <- function(gmt_path, prefix) {
  # 1. 读取基因集
  gs <- read.gmt(gmt_path)
  
  # 2. 基因集筛选（根据类型针对性处理）
  if (prefix %in% c("c2.all", "c2.kegg_legacy", "c2.kegg_medicus", "c2.reactome")) {
    # C2类：去除前缀（如"REACTOME_"），保留经典通路
    gs$term <- str_remove(gs$term, paste0("^", str_remove(prefix, "c2\\."), "_|^CP_"))
  } else if (prefix == "c5.all") {
    # C5（GO）：仅保留生物学过程（GOBP）
    gs <- gs %>% filter(str_detect(term, "GOBP_"))
    gs$term <- str_remove(gs$term, "GOBP_")
  } else if (prefix == "c7.all") {
    # C7（免疫）：筛选小胶质细胞/神经炎症相关通路
    keywords <- c("microglia", "neuroinflammation", "glial", "M1", "M2", "phagocyt")
    gs <- gs %>% filter(str_detect(tolower(term), paste(keywords, collapse = "|")))
  } else if (prefix == "h.all") {
    # Hallmark：去除前缀
    gs$term <- str_remove(gs$term, "HALLMARK_")
  }
  
  # 3. 转换为基因集列表并匹配表达矩阵基因
  gs_list <- gs %>% split(.$term) %>% lapply("[[", 2)
  gs_list <- lapply(gs_list, function(x) intersect(x, rownames(expr_matrix)))
  gs_list <- gs_list[sapply(gs_list, length) >= 2]  # 保留至少2个基因的通路
  if (length(gs_list) == 0) stop(paste(prefix, "无有效通路，检查基因匹配！"))
  
  # 4. 运行SSGSEA
  ssgsea_param <- ssgseaParam(
    exprData = expr_matrix,
    geneSets = gs_list,
    minSize = 2, maxSize = 500, alpha = 0.25, normalize = TRUE, verbose = FALSE
  )
  gsva_res <- gsva(ssgsea_param)
  
  # 5. 差异检验（Epilepsy vs Healthy）
  design <- model.matrix(~0 + factor(group_info))
  colnames(design) <- levels(factor(group_info))
  rownames(design) <- colnames(gsva_res)
  contrast <- makeContrasts(Epilepsy - Healthy, levels = design)
  fit <- lmFit(gsva_res, design) %>% contrasts.fit(contrast) %>% eBayes()
  diff_res <- topTable(fit, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  
  # 6. 绘制条形图
  plot_df <- data.frame(
    Pathway = rownames(diff_res),
    t_value = diff_res$t,
    FDR = diff_res$adj.P.Val
  ) %>% arrange(t_value) %>% mutate(Pathway = factor(Pathway, levels = Pathway))
  plot_df$signif <- ifelse(plot_df$FDR < 0.05 & abs(plot_df$t_value) > 2, "Significant", "Not")
  
  p <- ggplot(plot_df, aes(x = Pathway, y = t_value, fill = signif)) +
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("Significant" = "#C73E1D", "Not" = "#CCCCCC")) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", alpha = 0.7) +
    labs(x = "Pathways", y = "t-value (Epilepsy - Healthy)", fill = "Significant (FDR<0.05, |t|>2)") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom") +
    ggtitle(paste(prefix, "Pathway Differences"))
  
  # 7. 保存结果
  ggsave(paste0(prefix, "_SSGSEA_Barplot.png"), p, width = 10, height = 8, dpi = 300)
  write.csv(diff_res, paste0(prefix, "_SSGSEA_Diff.csv"), row.names = TRUE)
  mg[[paste0("GSVA_", prefix)]] <- CreateAssayObject(data = gsva_res)  # 添加到Seurat对象
  
  cat(paste0("\n✅ ", prefix, " 分析完成：有效通路数=", length(gs_list), "\n"))
  return(mg)  # 返回更新后的Seurat对象
}

# ------------------------------------------------------------------------------
# 2. 依次运行所有基因集分析（按gmt_files顺序）
for (i in seq_along(gmt_files)) {
  mg <- run_ssgsea(gmt_path = gmt_files[i], prefix = names(gmt_files)[i])
}

# 保存最终Seurat对象（含所有基因集的GSVA结果）
saveRDS(mg, "Microglia_with_All_GSVA.rds")

cat("\n所有基因集分析完成！结果文件：\n", 
    "- 每个基因集对应：条形图PNG + 差异通路CSV\n", 
    "- 合并Seurat对象：Microglia_with_All_GSVA.rds")