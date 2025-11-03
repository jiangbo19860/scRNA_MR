
BiocManager::install("GEOquery")
url='https://gitee.com/jmzeng/annoprobe.git'
library(remotes)
install_git(url)

library(Seurat)
library(infercnv)
load("scRNA_endo.rds")
DimPlot(scRNA)

scRNA$celltype
table(scRNA$celltype)

cells.use=colnames(scRNA)[which(scRNA$ClusterName %in% c("normal endothelial cell"))]

dat=as.data.frame(GetAssayData(subset(scRNA, cells=cells.use)))
groupinfo=data.frame(v1=colnames(dat),
                     v2=scRNA@active.ident[cells.use])

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match(geneInfor[,1], rownames(dat) ),] 
dim(dat)

groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)





infercnv_obj = CreateInfercnvObject(delim = '\t',
                                    raw_counts_matrix = 'expFile.txt',
                                    annotations_file = 'groupFiles.txt',
                                    gene_order_file = 'geneFile.txt',
                                    ref_group_names = c("Normal_KC_Cyc"))

#10x数据cutoff推荐使用0.1
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir='inferCNV/', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

##自己画矢量图
plot_cnv(infercnv_obj,
         out_dir="inferCNV2/",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         cluster_by_groups=TRUE,
         x.center=1,
         x.range="auto",
         hclust_method='ward.D',
         custom_color_pal = color.palette(c("#0071B2", "white", "#C3250A"), c(2, 2)),
         color_safe_pal=FALSE,
         output_filename="infercnv",
         output_format="pdf",
         #png_res=300,
         dynamic_resize=0)



