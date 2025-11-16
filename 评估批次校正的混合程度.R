# =============================================================================
# 评估批次校正的混合程度
# =============================================================================

# 1. 使用LISI评分（如果已安装lisi包）
if (require(lisi)) {
  # 计算每个细胞的样本混合度
  lisi_res <- compute_lisi(
    Embeddings(mg_filtered, "harmony")[, 1:30],
    mg_filtered@meta.data,
    c("Sample_ID")
  )
  
  # 添加到metadata
  mg_filtered$lisi_score <- lisi_res$Sample_ID
  
  # 按聚类可视化LISI分数
  VlnPlot(mg_filtered, 
          features = "lisi_score",
          group.by = "seurat_clusters",
          pt.size = 0) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "LISI Score by Cluster (Higher = Better Mixing)",
         caption = "Score close to 1 = poor mixing, close to # samples = good mixing")
  
  ggsave("lisi_by_cluster.pdf", width = 10, height = 6)
  
  # 识别混合不良的聚类
  lisi_by_cluster <- mg_filtered@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      mean_lisi = mean(lisi_score),
      median_lisi = median(lisi_score),
      n_cells = n()
    ) %>%
    arrange(mean_lisi)
  
  print(lisi_by_cluster)
  
  # 标记混合不良的聚类（LISI < 2表示主要来自单一样本）
  poorly_mixed <- lisi_by_cluster %>%
    filter(mean_lisi < 2) %>%
    pull(seurat_clusters)
  
  cat("\nPoorly mixed clusters (mean LISI < 2):\n")
  print(poorly_mixed)
}

# 2. 计算每个聚类的样本多样性（香农熵）
library(vegan)

calculate_diversity <- function(cluster_id) {
  cells_in_cluster <- mg_filtered@meta.data %>%
    filter(seurat_clusters == cluster_id)
  
  sample_counts <- table(cells_in_cluster$Sample_ID)
  shannon_entropy <- diversity(sample_counts, index = "shannon")
  
  return(shannon_entropy)
}

cluster_diversity <- sapply(levels(mg_filtered$seurat_clusters), 
                            calculate_diversity)

diversity_df <- data.frame(
  Cluster = levels(mg_filtered$seurat_clusters),
  Shannon_Entropy = cluster_diversity
) %>%
  arrange(Shannon_Entropy)

print(diversity_df)

# 可视化
ggplot(diversity_df, aes(x = reorder(Cluster, Shannon_Entropy), 
                         y = Shannon_Entropy)) +
  geom_col(aes(fill = Shannon_Entropy < 1.5)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Sample Diversity by Cluster (Shannon Entropy)",
       x = "Cluster",
       y = "Shannon Entropy",
       caption = "Lower entropy = dominated by fewer samples") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    name = "Low diversity\n(< 1.5)") +
  theme_classic()

ggsave("cluster_diversity.pdf", width = 10, height = 8)
