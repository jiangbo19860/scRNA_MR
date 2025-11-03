rm(list = ls())

# org.Hs.eg.db是人类基因组数据库，核心功能包括基因ID转换（Entrez ID ↔ 基因符号（Symbol）, Entrez ID ↔ Ensembl ID, Entrez ID ↔ UniProt ID）
pacman::p_load(here, Seurat, presto, fgsea, tidyverse, patchwork, clusterProfiler, AnnotationDbi, org.Hs.eg.db, fgsea, enrichplot, data.table, Matrix, ggplot2, ggrepel, ggpubr, pheatmap, ggsci, ggridges, cowplot, survminer, survival, limma, hdf5r)

here()

load(here("1_data", "E-MTAB-6149", "scRNA_endo.Rds"))

str(scRNA) # 查看scRNA对象的结构，重要的slots: @assays核心数据存储（包括counts, data, scale.data, var.features, meta.features），@meta.data细胞水平的元数据，是一个数据框，每行是一个细胞barcode,各种细胞的信息如nCount_RNA, nFeature_RNA, percent.mt等，@reductions降维和聚类结果（如PCA, tsne, UMAP等），@graphs图结构数据（如邻接矩阵，细胞临域图RNA_nn, 共享最近邻图RNA_snn, 是聚类分析的基础），@commands记录分析步骤（如标准化、找高变基因、PCA、聚类、降维等）。

dim(scRNA) # 查看scRNA对象的维度
DimPlot(scRNA, reduction = "umap", label = TRUE) # 绘制UMAP降维图

meta <- scRNA@meta.data


# 0. 差异基因比的是什么--把细胞分组---------------------------------------------------------------
# 查看默认的细胞身份标识（Identity）是什么及分布
table(Idents(scRNA))

# 为后续DEGs分析设置新的细胞分类标准：新增group列，按是否来自于肿瘤细胞分成两组，这些改变发生在meta.data中
scRNA$group <- ifelse(scRNA$CellFromTumor == "1", "Tumor", "Non-malignant") # 创建一个新的分组变量)

meta <- scRNA@meta.data

# 把细胞身份标识设置为上面设置的分组（按CellFromTumor分组）
Idents(scRNA) <- "group"
# 检查细胞身份标识是否设置成功
table(Idents(scRNA))
# 绘制UMAP降维图，检查细胞的分组是否按照预期分组
DimPlot(scRNA, reduction = "umap", label = TRUE)

# 1. 差异基因分析FindMarkers和火山图 -------------------------------------------------------------------

# FindMarkers函数识别差异表达基因,FindMarkers() 默认使用 Wilcoxon 秩和检验，presto 包提供了更高效的实现。
deg <- FindMarkers(scRNA,
  ident.1 = "Tumor",
  ident.2 = "Non-malignant"
)

colnames(deg) # 5列："p_val","avg_log2FC","pct.1","pct.2","p_val_adj".
# p_val：原始 p 值
# avg_log2FC：平均 log2 倍数变化（正值表示在 Tumor 中高表达，负值表示在 Non-malignant 中高表达）
# pct.1：在 Tumor 组中检测到该基因表达的细胞百分比
# pct.2：在 Non-malignant 组中检测到该基因表达的细胞百分比
# p_val_adj：校正后的 p 值（通常使用 Benjamini-Hochberg 方法）

# 量化在哪个组中是高表达，哪个组是低表达-用小提琴图VlnPlot可视化查看，形状：表示表达值的概率密度（越宽表示该区域细胞数量越多）。箱体：中间的箱体显示四分位数范围（25%-75%），中线为中位数。若小提琴图整体位置较高，说明该基因在多数细胞中表达较高；若贴近基线，则表达较低或不表达。若小提琴图形状宽且分散，说明细胞间表达差异大；若窄而集中，则差异小。若未指定 group.by 参数，Seurat 默认按当前身份标识（Idents(scRNA)）分组。

# 查看差异最大的一个基因在两组中的表达情况，结果显示在非肿瘤细胞中高表达，说明pt.2是非肿瘤细胞组，pt.1是肿瘤。avg_log2FC是肿瘤组细胞/非肿瘤组细胞。
VlnPlot(scRNA, features = "SLC25A6")

# 查看前10个差异基因在两组中的表达情况
VlnPlot(scRNA, features = rownames(deg)[1:10], group.by = "group") + NoLegend()

sig_deg <- subset(deg, p_val_adj < 0.01 & abs(avg_log2FC) > 1)

# 绘制差异基因火山图Volcano Plot，x轴表示基因表达的倍数变化（如 log2FC）, y轴表示校正后的 p 值（-log10(p_val_adj)）。
ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = factor(ifelse(p_val_adj < 0.01 & abs(avg_log2FC) > 1, "Significant", "Non-significant"),
    levels = c("Significant", "Non-significant")
  ))) +
  scale_color_manual(
    values = c("red", "gray"),
    labels = c("Significant", "Non-significant")
  ) +

  # 添加x轴和y轴的黑线（轴线）
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +

  # 添加阈值虚线（垂直虚线表示log2FC阈值，水平虚线表示p值阈值）
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted p-value)",
    color = "Group"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_classic() + # 使用经典主题，默认包含坐标轴
  theme(
    panel.background = element_blank(), # 背景透明
    axis.line = element_line(color = "black", linewidth = 0.5), # 确保坐标轴线条可见
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank() # 移除次要网格线
  )


# 2. GO富集分析---------------------------------------------------------------------
# 用clusterProfiler 包中的enrichGO()函数对差异表达基因进行 Gene Ontology (GO) 富集分析
ego_ALL <- enrichGO(
  gene = row.names(sig_deg), # 输入基因：差异表达基因的名称（行名）
  OrgDb = "org.Hs.eg.db", # 物种数据库：人类基因注释数据库
  keyType = "SYMBOL", # 基因ID类型：使用基因符号（Symbol）作为输入
  ont = "ALL", # ALL表示分析所有GO类别（BP, CC, MF），BP (Biological Process)：生物过程（如 “细胞凋亡”、“免疫应答”）。CC (Cellular Component)：细胞组分（如 “线粒体”、“细胞膜”）。MF (Molecular Function)：分子功能（如 “ATP 酶活性”、“DNA 结合”）。
  pAdjustMethod = "BH", # p 值校正方法：Benjamini-Hochberg 法控制 FDR
  pvalueCutoff = 0.01, # p 值阈值：0.01
  qvalueCutoff = 0.05 # q 值（BH校正后p值（adjusted p-value）阈值：0.05，控制假阳性率FDR的关键指标。
)

ego_all <- data.frame(ego_ALL) # 将富集结果转换为可操作的数据框
colnames(ego_all) # 有13列：
# 1. ONTOLOGY：GO术语的类别（BP生物过程，CC细胞组分，MF分子功能），
# 2. ID：GO 术语的唯一标识符（格式为 GO:XXXXXXX）。
# 3. Description：GO 术语的文字描述（如 “apoptotic process”），
# 4. GeneRatio：差异基因中属于该 GO 术语的基因比例。差异基因中注释到该 GO 术语的基因数/所有差异基因的总数。若 GeneRatio = 15/100，表示 100 个差异基因中有 15 个属于该 GO 术语。
# 5. BgRatio: 背景基因中属于该 GO 术语的基因比例。背景基因中注释到该 GO 术语的基因数/参考基因组的基因总数（通常为数据库收录的基因数）。若 BgRatio = 50/18888，表示18888个背景基因中有 50 个属于该 GO 术语。
# 6. RichFactor: 富集因子，计算公式为 RichFactor = GeneRatio / BgRatio。衡量差异基因在该 GO 术语中的富集程度，值越大表示富集越显著。若 RichFactor = 3，表示差异基因在该 GO 术语中的富集程度是基因组背景的 3 倍。
# 7. pvalue：富集分析的原始 p 值（通常基于超几何分布计算）。
# 8. p.adjust：Benjamini-Hochberg (BH) 校正后的 p 值，用于控制 FDR（假发现率）。解决多重检验中的假阳性问题。通常需 < 0.05，表示在所有被判定为显著的结果中，假阳性比例 ≤ 5%。
# 9. qvalue：clusterProfiler 中该列与 p.adjust 完全相同（实际是 BH 校正后的 p 值），并非严格意义上的 q 值。
# 10. FoldEnrichment: 富集倍数，计算公式为 FoldEnrichment = (GeneRatio / BgRatio) * 100。表示差异基因在该 GO 术语中的富集程度相对于背景基因的倍数。若 FoldEnrichment = 2，表示差异基因在该 GO 术语中的富集程度是背景基因的 2 倍。
# 11. zScore: z 分数，用于评估 GO 术语的表达偏离程度。z = (实际观察值 - 期望值) / 标准差。正值表示该 GO 术语在差异基因中高表达，负值表示低表达。
# 12. geneID: 属于该 GO 术语的差异基因列表（用分号;或/分隔）。如TP53;CASP3;BAX，表示这三个差异基因均注释到该 GO 术语。
# 13. Count：差异基因中属于该 GO 术语的基因数量（即 GeneRatio 中的分子 a）。若 Count = 15，表示 15 个差异基因注释到该 GO 术语。

#
ego_BP <- enrichGO(
  gene = row.names(sig_deg),
  # universe = row.names(dge.celltype),
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "BP", # 只分析生物过程(Biological Process)类别
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

# dotplot是clusterProfiler包中的气泡图函数，用于展示富集结果。
dotplot(ego_ALL, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")


# 3. KEGG富集分析 ----------------------------------------------------------------
# 基因ID类型转换（bitr: Biological ID TRanslator），因为KEGG仅识别Entrez ID，所以需要将差异表达基因的符号(SYMBOL)转换为Entrez ID。
genelist1 <- bitr(row.names(sig_deg),
  fromType = "SYMBOL", # 输入基因ID类型：基因符号（如"SLC25A6"）
  toType = "ENTREZID", # 输出基因ID类型：Entrez ID（如"6523"）
  OrgDb = "org.Hs.eg.db" # 人类基因注释数据库
)

genelist1 <- pull(genelist1, ENTREZID) # 从转换结果中提取Entrez ID列，得到一个基因ID向量。

# 执行KEGG富集分析
ekegg <- enrichKEGG(gene = genelist1, organism = "hsa") # hsa（Homo sapiens）表示人类基因组
# enrichKEGG函数用于对输入的基因进行KEGG通路富集分析，返回一个enrichResult对象。

dotplot(ekegg) # 绘制KEGG富集结果的点图，展示各通路的富集情况。
# 将富集结果转换为数据框，方便后续操作
kegg <- data.frame(ekegg)



# 4. GSEA (Gene Set Enrichment Analysis)-----------------------------------------------------------------
# 基于基因集的富集分析方法，通常用于分析差异表达基因集是否在预定义的基因集中富集。把所有差异基因按照平均log2FC从大到小排序，得到一个基因列表（geneList），其中基因名为Entrez ID，值为log2FC。

### GSEA第一步：创建geneList:包括3列：值：基因的avg_log2FC（对数倍变化，反映表达差异程度）。名称：对应的ENTREZID（基因唯一标识符）。顺序：按avg_log2FC从大到小排序（GSEA 需要按差异程度排序）。

# 进行ID转换(gene符号symbol转为ENTREZID)，保留所有基因（包括无法映射的，用NA标记），得到一个2列的数据框id_mapping（"SYMBOL"和"ENTREZID"）
id_mapping <- bitr(
  geneID = rownames(deg), # 输入原始基因名
  fromType = "SYMBOL", # 输入类型：基因符号
  toType = "ENTREZID", # 输出类型：Entrez ID
  OrgDb = org.Hs.eg.db,
  drop = FALSE # 关键：不删除无法映射的基因，用NA表示
)

# 原始数据deg新建ENTREZID列，将转换结果合并到原始数据框（保留行名与转换结果的对应关系）
deg$ENTREZID <- id_mapping$ENTREZID # 此时行数一致，未映射的为NA

# 删除无法映射的基因（ENTREZID为NA的行）
deg_filtered <- deg[!is.na(deg$ENTREZID), ] # !is.na()筛选出非NA的行

# 将基因名作为一列数据rownames_to_column(), 得到一个新的数据框deg_filtered，其中包含ENTREZID列和基因名列（gene）。
deg_filtered <- rownames_to_column(deg_filtered, var = "gene")

# 验证结果
nrow(deg) # 原始行数：10432
nrow(deg_filtered) # 过滤后行数：9695（与可映射的数量一致）
sum(is.na(deg_filtered$ENTREZID)) # 确认过滤后无NA：0

# 5. 为基因集富集分析（GSEA）准备输入数据：按差异表达程度(平均对数倍变化avg_log2FC)排序(从大到小)的基因列表（ENTREZID）。

# 先排序数据框
sorted_df <- deg_filtered %>%
  arrange(desc(avg_log2FC)) # 按avg_log2FC降序排列

# 从排序后的数据框中提取两列，创建命名向量
geneList2 <- setNames(sorted_df$avg_log2FC, sorted_df$ENTREZID)

# 检查向量结构
str(geneList2) # 应显示 "Named num [1:?]"
head(geneList2)

# 可视化命名向量
# 将命名向量转换为数据框
df <- data.frame(
  ENTREZID = names(geneList2),
  avg_log2FC = as.numeric(geneList2)
)

# 绘制柱状图，x轴为基因ID，y轴为log2FC
ggplot(df, aes(x = ENTREZID, y = avg_log2FC)) +
  geom_col() +
  labs(title = "Gene Expression Changes", x = "ENTREZ ID", y = "avg_log2FC")

gsea <- gseKEGG(geneList = geneList2, organism = "hsa") # 执行GSEA分析，使用KEGG数据库
dotplot(gsea)

# 使用更美观的可视化
dotplot(gsea, showCategory = 20) +
  theme_minimal() +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    title = "GSEA Enrichment of KEGG Pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "KEGG Pathway"
  )

# 查看显著富集的通路数量
summary(gsea)

# 提取显著通路并排序
sig_pathways <- gsea[gsea$p.adjust < 0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES, decreasing = TRUE), ]

# 查看富集程度最高的通路
head(sig_pathways)


# 自定义基因集 ------------------------------------------------------------------
# 获得基因集：https://www.gsea-msigdb.org/gsea/msigdb/index.jsp -- Molecular Signatures Database -- Human MSigDB Collections -- 下载文件，如h.all.v2025.1.Hs.symbols.gmt文件，包含人类基因集的HALLMARK通路信息。-- 选出自己想要的，写成一个txt文件，然后保存为.gmt格式。

### 使用fgsea 包对差异表达基因进行 基因集富集分析（GSEA） ------------------------------------------

# 准备排序后的基因表达值
# 更简洁的写法
head(rownames(deg)) # 应显示基因符号（如"SFTPC"    "MT2A"     "HSPG2"    "PLVAP"    "ZFP36"    "HLA-DRB1"）
ranked_log2FC <- deg %>%
  arrange(desc(avg_log2FC)) %>%
  {
    setNames(.$avg_log2FC, rownames(.))
  } %>% # 使用行名（基因符号）
  .[!is.na(names(.))] %>% # 过滤缺失名称的基因
  .[!duplicated(names(.))] # 过滤重复名称的基因

str(ranked_log2FC) # 应显示 "Named num [1:?]"
head(ranked_log2FC) # 应显示基因名称和对应值

hallmark <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt")) # 这里的symbols要与上面的setNames中作为向量名称的基因符号一致
colnames(hallmark)
head(hallmark)
hallmark$term <- gsub("HALLMARK_", "", hallmark$term)

# 转换为列表格式（pathway名称 → 基因列表）
hallmark.list <- hallmark %>%
  split(.$term) %>%
  lapply("[[", 2)

# 验证hallmark.list是否正确创建
head(hallmark.list$APOPTOSIS) # 应显示凋亡通路的基因

# 执行GSEA分析，得到fgseaRes数据框, 有8列，[1] "pathway","pval","padj", "log2err","ES", "NES", "size","leadingEdge".
fgseaRes <- fgsea(
  pathways = hallmark.list,
  stats = ranked_log2FC,
  minSize = 15,
  maxSize = 500
)

colnames(fgseaRes)

# 按p值排序并查看前10个通路
sig_pathways <- fgseaRes[order(fgseaRes$padj), ]
head(sig_pathways, 10)

# 筛选显著通路（FDR < 0.25）
sig_pathways <- sig_pathways[sig_pathways$padj < 0.25, ]
nrow(sig_pathways) # 显著通路数量

# 例如，查看凋亡通路（HALLMARK_APOPTOSIS）
apoptosis <- fgseaRes[fgseaRes$pathway == "APOPTOSIS", ]
print(apoptosis)
class(apoptosis)

# 查看该通路的Leading Edge基因
leading_genes <- apoptosis$leadingEdge[[1]] # 提取列表的第一个元素
print(leading_genes) # 对凋亡通路富集贡献最大的基因

deg[leading_genes, "avg_log2FC"] # 假设deg是差异表达数据框

# 可视化单个基因集在差异基因中的富集情况曲线图，展示基因集（如凋亡通路）的基因在排序基因列表中的分布模式。绘制富集曲线并标记Leading Edge， 横坐标（X 轴）：rank(按 log2FC（差异表达倍数的对数值） 或 t值 排序), 代表所有的deg基因按差异程度排序后的位置（从左到右，基因差异程度逐渐变化）。左侧（低 rank）：显著上调基因（高表达，与表型正相关）。右侧（高 rank）：显著下调基因（低表达，与表型负相关）。纵坐标（Y 轴）：enrichment score, 代表 GSEA 算法滑动窗口计算的富集得分（反映通路基因在排序基因中的 “聚集程度”）。曲线先上升后下降（或相反），峰值位置对应通路基因最富集的区域（即对通路贡献最大的基因子集，Leading Edge）。绿色曲线（Enrichment Score 轨迹）
plotEnrichment(hallmark.list[["APOPTOSIS"]], ranked_log2FC) +
  labs(title = "Apoptosis Pathway Enrichment")


# 添加Leading Edge基因标记
# 示例：只标记NES贡献最大的前20个基因
# 1. 计算Leading Edge基因的log2FC绝对值（反映重要性）
leading_log2fc <- ranked_log2FC[leading_genes]
gene_importance <- abs(leading_log2fc)

# 2. 筛选前20个最重要的基因
top_genes <- names(sort(gene_importance, decreasing = TRUE)[1:20])
top_positions <- match(top_genes, names(ranked_log2FC))

# 3. 绘图
p +
  geom_vline(xintercept = top_positions, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text_repel(
    data = data.frame(x = top_positions, y = 0.5, label = top_genes),
    aes(x = x, y = y, label = label),
    color = "red",
    size = 3,
    max.overlaps = 20 # 与筛选数量一致
  )

# quartzFonts()
# sig <- fgseaRes
# sig <- sig[order(sig$NES, decreasing = T), ]
# up <- sig[1:10]$pathway
# sig <- sig[order(sig$NES, decreasing = F), ]
# down <- sig[1:10]$pathway
# ## 最后一步 开始绘图
# library(ggplot2)
# library(fgsea)
#
# # 全局设置通用无衬线字体（避免直接指定 Helvetica）
# theme_set(theme_minimal(base_family = "sans"))
#
# # 重新绘图
# p2 <- plotGseaTable(
#   hallmark.list[c(up, down)],
#   id,
#   fgseaRes,
#   gseaParam = 0.5
# )
# p2


# 用clusterProfiler::GSEA()进行分析，需提供TERM2GENE格式，适合自定义基因集 --------------------
# 对所有差异表达基因（deg）按平均对数倍变化（avg_log2FC）降序排列。
alldiff <- deg[order(deg$avg_log2FC, decreasing = T), ]

# 创建 GSEA 所需的命名数值向量（基因名→log2FC 值）。

# 提取排序后的 log2FC 列。
id <- alldiff$avg_log2FC

# 为 log2FC 值添加基因名作为名称。id 是一个命名向量，如 c("EGFR"=3.2, "KRAS"=2.8, ...)。
names(id) <- rownames(alldiff)
# 如果是使用 ENTREZID 列作为基因名
# names(id) <- alldiff$ENTREZID  # 假设 alldiff 包含 ENTREZID 列

head(id)

hallmark <- read.gmt(here("1_data", "GSEA", "h.all.v2025.1.Hs.symbols.gmt"))

# 执行 GSEA 分析，使用 fgsea 包, id是命名的 log2FC 向量，TERM2GENE：基因集映射关系。
y <- GSEA(id, TERM2GENE = hallmark)

# 绘制气泡图展示富集结果,横轴为通路名，纵轴为富集分数，气泡大小表示基因数，颜色表示 FDR。showCategory = 12：显示前 12 个通路。split = ".sign"：按富集方向（上调 / 下调）分组。facet_grid(~.sign)：按富集方向分面展示。
dotplot(y, showCategory = 12, split = ".sign") + facet_grid(~.sign)

# 将 GSEA 结果转换为可操作的数据框。
yd <- data.frame(y)

# 可视化某个具体通路（如HALLMARK_HYPOXIA） 的富集情况。y：GSEA 结果对象。"HALLMARK_HYPOXIA"：指定要展示的通路名称。color = "red"：曲线颜色为红色。pvalue_table = T：在图下方显示 p 值和 NES 表格。
gseaplot2(y, "HALLMARK_HYPOXIA", color = "red", pvalue_table = T)
