rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggpubr)

load("Fibro.Cellview.Rds")
log2cpm[1:3,1:3]
load("scRNA.Rds")

scRNA<-scRNA[,colnames(log2cpm)]
scRNA <-AddMetaData(scRNA,metadata = tsne.data)
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000) 
##选择高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# top10 <- head(VariableFeatures(scRNA), 10) 


##数据中心化与回归
#数据中心化
#scale.genes <-  rownames(scRNA)
#scale.genes <-  VariableFeatures(scRNA)    #内存不够可用此命令代替上一行
scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 

DimPlot(scRNA, reduction = "pca", group.by="orig.ident")

ElbowPlot(scRNA, ndims=20, reduction="pca") 

pc.num=1:15

scRNA <- FindNeighbors(scRNA, dims = pc.num) 

scRNA <- FindClusters(scRNA, resolution = 0.2)

metadata <- scRNA@meta.data

scRNA = RunTSNE(scRNA, dims = pc.num)
scRNA <- RunUMAP(scRNA, dims = pc.num)

save(scRNA,file = "scRNA_fibro.Rds")

load("scRNA_fibro.Rds")

DimPlot(scRNA,group.by = "dbCluster",reduction = "tsne")

DimPlot(scRNA,reduction = "tsne")

DimPlot(scRNA,reduction = "tsne",label = T)+
  theme(panel.border = element_rect(fill=NA,color="black", 
                                    size=1, linetype="solid"),
        legend.position = "none")+
  labs(title = "Fibroblasts")

ggsave(filename = "Fig/Fig3a.pdf",height = 8,width = 8)

####################################################################################

DimPlot(scRNA, 
        reduction = "tsne",
        group.by ="CellFromTumor")+
  scale_color_manual(values = c("#72C667","#2D5474"),
                     labels=c("Non-malignant","Tumor"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")

ggsave(filename = "Fig/Fig3a1.pdf",height = 8,width = 8)


####################################################################################
tmp=scRNA@meta.data
tmp$group=ifelse(tmp$CellFromTumor=="1","Tumor","Non-malignant")
tmp$seurat_clusters =factor(tmp$seurat_clusters,levels = 5:0)

p1=ggplot(tmp,
          aes(x=seurat_clusters,
              fill=group))+
  geom_bar(stat="count",position="fill")+
  theme_classic()+
  scale_fill_manual(values=c("#72C667","#2D5474"))+
  theme(axis.text.x = element_text(angle=90,
                                   vjust=0.5))+xlab(" ")+ylab(" ")+
  theme(legend.position = "top")+
  theme(axis.line.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.y = element_blank())+
  coord_flip()

p2=ggplot(tmp,
          aes(x=seurat_clusters,fill=PatientNumber))+
  geom_bar(stat="count",position="fill")+
  theme_classic()+
  scale_fill_manual(values=c("#268A24",
                             "#8EEE8B",
                             "#FB6346",
                             "#FBD51A",
                             "#28507D"))+
  theme(axis.text.x = element_text(angle=90,
                                   vjust=0.5))+xlab(" ")+ylab(" ")+
  theme(legend.position = "top")+
  theme(axis.line.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.y = element_blank())+
  coord_flip()
p1+p2+plot_layout(ncol = 2)

ggsave("Fig/Fig3_fibro_bar.pdf",height = 4,width = 10)


################################################################################

p=FeaturePlot(scRNA, 
              reduction = "tsne",
              features = c("COL10A1",
                           "COL4A1",
                           "PLA2G2A",
                           "MMP3",
                           "FIGF",
                           "CCL2"),
              ncol = 3,
              cols = c("gray","red"))
p & theme(panel.border = element_rect(fill=NA,color="black", 
                                      size=1, linetype="solid"),
          legend.position = "none")

ggsave("Fig/Fig3b.pdf",height = 6,width = 9)

######################################################################################

library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

gs <-read.gmt("h.all.v7.4.symbols.gmt")

gs$term <- gsub('HALLMARK_','',gs$term)

gs.list <- gs %>% split(.$term) %>% lapply( "[[", 2)

gsva_es <- gsva(as.matrix(scRNA@assays$RNA@counts), 
                gs.list,
                method="ssgsea",
                abs.ranking=T)

data_gs <-data.frame(row.names = rownames(gsva_es))

for (x in 0:5) {
  group <-rownames(scRNA@meta.data)[which(scRNA$seurat_clusters==x)]
  gsva_es_tmp <-gsva_es[,group]
  gsva_es_tmp <-apply(gsva_es_tmp, 1, mean)
  gsva_es_tmp <-as.data.frame(gsva_es_tmp)
  colnames(gsva_es_tmp)<-x
  data_gs <-cbind(data_gs,gsva_es_tmp)
}

pdf("Fig/Fig3f.pdf",height = 10,width = 8)
library(pheatmap)
pheatmap(data_gs,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = F,
         cluster_rows = F,
         scale = "row")
dev.off()





data=data.frame(row.names = rownames(gsva_es))

for (i in 0:5) {
  scRNA$group <-ifelse(scRNA$seurat_clusters==i,
                       "A",
                       "B")
  
  group_list <- data.frame(sample = colnames(gsva_es), 
                           group = scRNA$group)
  
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  contrast.matrix <- makeContrasts(A-B, levels = design)
  fit <- lmFit(gsva_es, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  data_x=data.frame(row.names = rownames(gsva_es),
                    t=as.numeric(x$t))
  colnames(data_x)=i
  data=cbind(data,data_x)
}

library(pheatmap)
pheatmap(data,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = F,
         cluster_rows = F,
         scale = "row")


#################################################################################


gene_name<-c("COL6A2",
             "COL4A1",
             "COL14A1",
             "COL5A1",
             "COL5A2",
             "COL8A1")
source("stackvlion.R")
StackedVlnPlot(scRNA,features = gene_name,idents=c(5,1,4,3,0)) & ggsci::scale_fill_lancet()

ggsave(filename = "Fig/Fig3c.pdf",height = 6,width = 5)












