rm(list = ls())
# remotes::install_github("WangJin93/GCAS")
library(GCAS)
get_expr_data #  get_expr_data() 是 GCAS 包中用于批量获取基因表达量数据的工具函数。get_expr_data() 函数的源代码、字节码位置和所属命名空间，核心功能是从指定数据集中获取基因的表达量数据，使用时只需传入数据集名称（datasets）和基因列表（genes），即可得到包含表达量和样本信息的整合结果，方便后续进行差异表达分析、可视化等操作。

results <- get_expr_data(datasets = "GSE74706", genes = c("GAPDH","TNS1"))
single_gene_multiple_datasets <- get_expr_data(datasets = c("GSE62113","GSE74706"), genes = "GAPDH")
multiple_genes_multiple_datasets <- get_expr_data(datasets = c("GSE62113","GSE74706"), genes = c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA"))

# 可视化GEO数据库中肿瘤组织和正常组织之间的mRNA表达数据的不同。
df_single <- get_expr_data(datasets = "GSE27262",genes = c("TP53"))
viz_TvsN(df_single,df_type = "single")

df_multi_gene <- get_expr_data(datasets = "GSE27262",genes = c("TP53","TNS1"))
viz_TvsN(df_multi_gene,df_type = "multi_gene",tumor_subtype ="LC")

df_multi_set <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
viz_TvsN(df_multi_set,df_type = "multi_set")

# 计算不同数据集中基因表达数据的摘要统计量（均值、标准差等）并进行假设检验（t 检验或 Wilcoxon 检验）。
df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")

# 绘制CPTAC数据集中肿瘤样本与正常样本之间差异表达基因（DEGs）的火山图。该功能对多个数据集进行荟萃分析，并生成森林图。同时，它还测试出版偏倚。
df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
plot_meta_forest(results)

# 生成基因在不同数据集中的对数倍数变化（log fold change, logFC）热图。热图中包含基于p值的显著性注释。
df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
heatmap <- plot_logFC_heatmap(results)
print(heatmap)

# 生成基因在不同数据集中的对数倍数变化（log fold change, logFC）散点图。散点图中包含基于p值的显著性注释。
df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
scatter <- plot_logFC_scatter(results, logFC.cut = 0.5, colors = c("blue","grey20", "red"))
print(scatter)

# 对CPTAC数据库中的mRNA/蛋白质表达数据执行相关性分析。
results <- cor_cancer_genelist(dataset = "GSE62113",
                               id1 = "STAT3",tumor_subtype = "LC",
                               id2 = c("TNS1", "TP53"),
                               sample_type = c("Tumor", "Normal"),
                               cor_method = "pearson")

# 计算目标基因表达与抗肿瘤药物敏感性在多个数据集之间的相关性。
dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210",
             "GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072",
             "GSE74706","GSE18842","GSE62113")
df <- get_expr_data(genes = "TNS1", datasets = dataset)
result <- cor_gcas_drug(df, Target.pathway = c("Cell cycle"))

# 对多个数据集中的表达数据执行相关性分析。
genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))

# 计算目标基因表达与免疫细胞浸润在多个数据集之间的相关性。
dataset <- c("GSE27262", "GSE7670", "GSE19188", "GSE19804", "GSE30219",
             "GSE31210", "GSE32665", "GSE32863", "GSE43458", "GSE46539",
             "GSE75037", "GSE10072", "GSE74706", "GSE18842", "GSE62113")
df <- get_expr_data(genes = "TNS1", datasets = dataset)
result <- cor_gcas_TIL(df, cor_method = "spearman", TIL_type = "TIMER")

# 使用基于 ggplot2 的热图展示相关性分析结果。
viz_cor_heatmap(result$r, result$p)

# 绘制散点图，包含样本大小（n）、相关系数（r）和p值（p.Value）。
viz_corplot(result$sss$GSE10072,"T_cell_CD4_TIMER","TNS1",x_lab = "")

