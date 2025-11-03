library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
rm(list=ls())

load("scRNA_endo.Rds")
DimPlot(scRNA)


scRNAsub <-subset(scRNA,idents=c(0,1,2))

DimPlot(scRNAsub)

##创建monocle的CDS对象
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        #lowerDetectionLimit = 1,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##选择用于降维排序的基因（三选一）
# 使用clusters差异表达基因######################################################

# diff.genes <- read.csv('subcluster_diff_genes_wilcox.csv')
# diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
# mycds <- setOrderingFilter(mycds, diff.genes)
# plot_ordering_genes(mycds)

#使用seurat选择的高变基因#######################################################

# var.genes <- VariableFeatures(scRNA)
# mycds <- setOrderingFilter(mycds, var.genes)
# plot_ordering_genes(mycds)


# 使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)


##降维及拟时分析
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot2
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plot3
##合并作图
plot1|plot2|plot3
pheno <-Biobase::pData(mycds)

scRNAsub$pseudotime <-pheno$Pseudotime

FeaturePlot(scRNAsub,features = "pseudotime")

plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 1)
plot_cell_trajectory(mycds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 1)

##monocle基因可视化
s.genes <- c("EPCAM")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
plotc

##拟时相关基因聚类热图
#cluster差异基因
diff.genes <- FindAllMarkers(scRNAsub,only.pos = T)

sig_diff.genes <- subset(diff.genes,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene

sig_diff.genes <- unique(as.character(sig_diff.genes))

diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

sig_gene_names <- row.names(subset(diff_test, qval < 0.01))

plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=3,
                             show_rownames=T, return_heatmap=T)

#高变基因
disp_table <- dispersionTable(mycds)

disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)

disp.genes <- as.character(disp.genes$gene_id)

diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))

plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)


##细胞轨迹分支相关基因
disp_table <- dispersionTable(mycds)

disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)

disp.genes <- as.character(disp.genes$gene_id)

mycds_sub <- mycds[disp.genes,]

plot_cell_trajectory(mycds_sub, color_by = "State")

beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)

beam_res <- beam_res[order(beam_res$qval),]

beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]

mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]

plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T)

##保存结果
saveRDS(mycds,file = "mycds.rds")







