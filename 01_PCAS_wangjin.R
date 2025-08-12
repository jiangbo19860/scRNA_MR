rm(list = ls())
remotes::install_github("WangJin93/PCAS")
# 加载PCAS包（需先安装）
library(PCAS)

# 1. get_data(): 获取CPTAC数据（通过API连接MySQL数据库）
# 示例：获取LUAD数据集的蛋白表达数据（GAPDH和TNS1基因）
expr_data <- get_data(
  table = "LUAD_Academia_protein",  # 数据集表格名称
  action = "expression",            # 操作类型：表达量数据
  genes = c("GAPDH", "TNS1")        # 目标基因符号
)
# 查看结果
head(expr_data)

# 2. get_expr_data(): 获取CPTAC数据库中的mRNA/蛋白表达数据
# 示例：从多个数据集获取TP53和TNS1的表达量
multi_expr <- get_expr_data(
  datasets = c("LUAD_CPTAC_protein", "LSCC_CPTAC_protein"),  # 数据集名称
  genes = c("TP53", "TNS1")                                 # 目标基因
)
# 查看结果
head(multi_expr)

# 3. Get_DEGs_result(): 获取肿瘤与正常样本的差异表达基因结果
# 示例：对LUAD蛋白数据集进行t检验分析差异基因
degs_result <- get_DEGs_result(
  dataset = "LUAD_CPTAC_protein",  # 数据集
  method = "t.test"                # 分析方法（t检验或limma）
)
# 查看结果
head(degs_result)

# 4. merge_clinic_data(): 合并表达数据与临床数据
# 示例：将表达数据与LUAD队列的临床数据合并
# 先获取表达数据
expr_data <- get_expr_data(datasets = "LUAD_CPTAC_protein", genes = c("TP53"))
# 合并临床数据
expr_clinic <- merge_clinic_data(
  cohort = "LUAD_CPTAC",  # 队列名称
  data_input = expr_data  # 表达数据（来自get_expr_data的输出）
)
# 查看结果
head(expr_clinic)

# 5. cor_cancer_genelist(): 单个数据集中基因表达相关性分析
# 示例：分析LUAD蛋白数据中STAT3与TNS1、TP53的相关性
cor_result <- cor_cancer_genelist(
  dataset1 = "LUAD_CPTAC_protein",  # 第一个数据集
  id1 = "STAT3",                    # 基因1
  dataset2 = "LUAD_CPTAC_mRNA",     # 第二个数据集
  id2 = c("TNS1", "TP53"),          # 基因2（可多个）
  sample_type = c("Tumor", "Normal"),  # 样本类型
  cor_method = "pearson"              # 相关性方法
)
# 查看结果
print(cor_result)

# 6. cor_pancancer_genelist(): 多数据集中基因表达相关性分析
# 查看PCAS包支持的所有数据集（通常通过内置的dataset对象获取）
?get_expr_data
results <- get_expr_data(datasets = c("CCRCC_CPTAC_mRNA","GBM_CPTAC_mRNA","HNSCC_CPTAC_mRNA","LSCC_CPTAC_mRNA","LUAD_CPTAC_mRNA","PDAC_CPTAC_mRNA","UCEC_CPTAC2_mRNA","UCEC_CPTAC1_mRNA"), genes = c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA"))
class(results)
head(results)

# 1. 获取多数据集的基因表达数据
gene1_expr <- get_expr_data(datasets = c("LUAD_CPTAC_mRNA", "LSCC_CPTAC_mRNA"), genes = "TP53")  # 目标基因（单个）
gene2_expr <- get_expr_data(datasets = c("LUAD_CPTAC_mRNA", "LSCC_CPTAC_mRNA"), genes = c("TNS1", "GAPDH"))  # 基因列表（多个）

# 2. 整合表达数据（保留此步骤，确保样本匹配）
df <- merge(gene1_expr, gene2_expr, by = c("ID", "dataset", "type"), all = TRUE)

# 3. 无需定义基因名称列表，直接使用gene2_expr作为geneset_data
# 4. 正确调用相关性分析函数（修正geneset_data参数）
cor_pancan_result <- cor_pancancer_genelist(
  df = df,                  # 整合后的表达数据（目标基因+基因列表）
  geneset_data = gene2_expr,  # 基因列表的表达数据（直接使用get_expr_data的输出）
  sample_type = c("Tumor", "Normal"),
  cor_method = "spearman"
)

# 查看结果
head(cor_pancan_result)

# 7. cor_pancancer_drug(): 基因表达与药物敏感性相关性分析
dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein","UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
result <- cor_pancancer_drug(df,Target.pathway = c("Cell cycle"))
# 查看结果
print(drug_cor)
head(drug_info)

# 8. cor_pancancer_TIL(): 基因表达与免疫细胞浸润相关性分析
dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein","LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein","UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
result <- cor_pancancer_TIL(df, cor_method = "spearman", TIL_type = c("TIMER"))


# 9. viz_TvsN(): 可视化肿瘤与正常样本的表达差异
# 示例：可视化TP53在肿瘤和正常样本中的表达差异
# 先获取数据
tp53_expr <- get_expr_data(datasets = "LUAD_CPTAC_protein", genes = "TP53")
# 绘图
viz_TvsN(
  df = tp53_expr,               # 表达数据
  df_type = "single",           # 单基因模式
  Show.P.value = TRUE,          # 显示P值
  Method = "t.test",            # 统计方法
  values = c("#00AFBB", "#FC4E07")  # 颜色设置
)

# 10. viz_DEGs_volcano(): 绘制差异基因火山图
# 示例：对差异分析结果绘制火山图
# 使用之前的差异分析结果degs_result
viz_DEGs_volcano(
  df = degs_result,             # 差异分析结果
  p.cut = 0.05,                 # P值阈值
  logFC.cut = 1,                #  Fold Change阈值
  show.labels = c("TP53", "TNS1")  # 标注特定基因
)

# 11. viz_cor_heatmap(): 相关性热图可视化
genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
dataset <- c("CCRCC_CPTAC_mRNA","GBM_CPTAC_mRNA","HNSCC_CPTAC_mRNA","LSCC_CPTAC_mRNA","LUAD_CPTAC_mRNA","PDAC_CPTAC_mRNA","UCEC_CPTAC2_mRNA","UCEC_CPTAC1_mRNA")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
result <- cor_pancancer_genelist(df,geneset_data,sample_type = c("Tumor"))
viz_cor_heatmap(result$r,result$p)


# 12. viz_corplot(): 绘制相关性散点图（带统计信息）
# 示例：绘制TP53与TNS1的表达相关性散点图
# 从multi_expr中提取两基因的表达数据
cor_data <- multi_expr[, c("TP53", "TNS1", "ID", "dataset")]
viz_corplot(
  data = cor_data,              # 包含两个基因表达数据的数据框
  a = "TP53",                   # 基因A
  b = "TNS1",                   # 基因B
  method = "pearson",           # 相关性方法
  x_lab = "TP53 expression",    # X轴标签
  y_lab = "TNS1 expression"     # Y轴标签
)

# 13. viz_phoso_sites(): 查询并可视化蛋白质磷酸化位点
# 示例：查询YTHDC2的磷酸化位点（基于CPTAC数据库）
viz_phoso_sites(
  gene = "YTHDC2",              # 目标蛋白
  phoso_infoDB = "CPTAC"        # 数据库选择（CPTAC或UniProt）
)

