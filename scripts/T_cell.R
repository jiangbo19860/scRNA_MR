rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)

load("T_cell.Cellview.Rds")
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

scRNA <- RunTSNE(scRNA, dims = pc.num)

scRNA <- RunUMAP(scRNA, dims = pc.num)

save(scRNA,file = "scRNA_T_cell.Rds")

load("scRNA_T_cell.Rds")

DimPlot(scRNA,group.by = "dbCluster",reduction = "tsne")

DimPlot(scRNA,reduction = "tsne")

DimPlot(scRNA,reduction = "tsne",label = T)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "none")+
  labs(title = "T cell")

ggsave(filename = "Fig/Fig5a.pdf",height = 8,width = 8)

####################################################################################

DimPlot(scRNA, 
        reduction = "tsne",
        group.by ="CellFromTumor")+
  scale_color_manual(values = c("#72C667","#2D5474"),
                     labels=c("Non-malignant","Tumor"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")

ggsave(filename = "Fig/Fig5a1.pdf",height = 8,width = 8)


####################################################################################
tmp=scRNA@meta.data
tmp$group=ifelse(tmp$CellFromTumor=="1","Tumor","Non-malignant")
tmp$seurat_clusters =factor(tmp$seurat_clusters,levels = 7:0)

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

ggsave("Fig/Fig5_T_cell_bar.pdf",height = 4,width = 10)


################################################################################

p=FeaturePlot(scRNA, 
              reduction = "tsne",
              features = c("FGFBP2",
                           "CD8A",
                           "CD4",
                           "FOXP3",
                           "LAG3",
                           "MKI67"),
              ncol = 3,
              cols = c("gray","red"))
p & theme(panel.border = element_rect(fill=NA,color="black", 
                                      size=1, linetype="solid"),
          legend.position = "none")

ggsave("Fig/Fig5b.pdf",height = 6,width = 9)

######################################################################################
Idents(scRNA)<-"ClusterName"
DimPlot(scRNA,reduction = "tsne",label = T)
library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

gs <-read.gmt("h.all.v7.4.symbols.gmt")

gs$term <- gsub('HALLMARK_','',gs$term)

gs.list <- gs %>% split(.$term) %>% lapply( "[[", 2)

clusternames=unique(scRNA@active.ident)

data_tmp<-data.frame(row.names = unique(gs$term))

for (i in 1:4) {
  scRNAsub <-subset(scRNA, ClusterName==clusternames[i])

  gsva_es <- gsva(as.matrix(scRNAsub@assays$RNA@counts), 
                  gs.list,
                  method="ssgsea",
                  abs.ranking=T)
  
  scRNAsub$group <-ifelse(scRNA$CellFromTumor=="1","tumor","Non_malignant")
  
  group_list <- data.frame(sample = colnames(gsva_es), 
                           group = scRNAsub$group)
  
  head(group_list)
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  design
  contrast.matrix <- makeContrasts(tumor-Non_malignant, levels = design)
  fit <- lmFit(gsva_es, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = 1, n = Inf, 
                adjust.method = "BH", 
                sort.by = "P")
  x <- data.frame(row.names = rownames(x),
                  t=x$t)
  colnames(x)<- clusternames[i]
  data_tmp<-cbind(data_tmp,x)
}


write.csv(data_tmp, "T_cell_gsva_limma.csv", quote = F)

library(pheatmap)
pdf("Fig/Fig5d.pdf",height = 10,width = 5)
pheatmap(data_tmp,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = F,
         cluster_rows = F
         )
dev.off()
#################################################################################

## CD8+ T cells

#################################################################################

meta=scRNA@meta.data
unique(meta$ClusterName)

scRNA=subset(scRNA,ClusterName=="CD8+ T cells ")

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

scRNA <- FindClusters(scRNA, resolution = 0.1)

metadata <- scRNA@meta.data

scRNA <- RunTSNE(scRNA, dims = pc.num)

scRNA <- RunUMAP(scRNA, dims = pc.num)

save(scRNA,file = "scRNA_CD8T_cell.Rds")

DimPlot(scRNA,group.by = "dbCluster",reduction = "tsne")

DimPlot(scRNA,reduction = "tsne")

library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

gs <-read.gmt("h.all.v7.4.symbols.gmt")

gs$term <- gsub('HALLMARK_','',gs$term)

gs.list <- gs %>% split(.$term) %>% lapply( "[[", 2)

data_tmp<-data.frame(row.names = unique(gs$term))

for (i in 0:3) {
  scRNAsub <-subset(scRNA,seurat_clusters==i)
  
  gsva_es <- gsva(as.matrix(scRNAsub@assays$RNA@counts), 
                  gs.list,
                  method="ssgsea",
                  abs.ranking=T)
  
  scRNAsub$group <-ifelse(scRNA$CellFromTumor=="1","tumor","Non_malignant")
  
  group_list <- data.frame(sample = colnames(gsva_es), 
                           group = scRNAsub$group)
  
  head(group_list)
  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_es)
  design
  contrast.matrix <- makeContrasts(tumor-Non_malignant, levels = design)
  fit <- lmFit(gsva_es, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = 1, n = Inf, 
                adjust.method = "BH", 
                sort.by = "P")
  x <- data.frame(row.names = rownames(x),
                  t=x$t)
  colnames(x)<- paste0("CLuster",i)
  data_tmp<-cbind(data_tmp,x)
}

library(pheatmap)
pdf("Fig/Fig5e.pdf",height = 10,width = 5)
pheatmap(data_tmp,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = F,
         cluster_rows = F,
         scale = "row")
dev.off()

################################################################################


gene_name<-c("PDCD1",
             "CD28",
             "CTLA4",
             "ICOS",
             "BTLA",
             "LAG3",
             "TNFRSF9",
             "TNFRSF4",
             "CD27",
             "HAVCR2",
             "GZMA",
             "GZMB",
             "GZMM",
             "GZMH",
             "GZMK")
source("stackvlion.R")
StackedVlnPlot(scRNA,features = gene_name,idents=c(0,1,2,3)) & ggsci::scale_fill_lancet()

ggsave(filename = "Fig/Fig5e2.pdf",height =15,width =4)


################################################################################
FeaturePlot(scRNA,features = c("GZMA",
                               "GZMB",
                               "GZMH"))



exp<-GetAssayData(scRNA)
exp<-as.data.frame(exp)
exp<-exp[,colSums(exp)>0]
exp<-exp[rowSums(exp)>0,]

GZ <-exp[c("GZMA","GZMB","GZMH"),]
GZ <-as.data.frame(GZ)
for (i in 1:ncol(GZ)) {
  GZ["GZ",i]<-mean(as.numeric(GZ[1:3,i]))
}

y <- as.numeric(GZ[4,])
rownames <- rownames(exp)
cor <-do.call(rbind,lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exp[x,]),y,type="spearman")
    data.frame(gene=x,
               cor=dd$estimate,
               p.value=dd$p.value)
  }))

save(cor,file = "T_cor.Rds")
cor=na.omit(cor)
cor=cor[order(cor$cor,decreasing = T),]
top=30
cor$label=c(cor$gene[1:top],rep(NA,(nrow(cor)-top)))
# cor$p=ifelse(cor$cor>=0,2-(cor$p.value),(cor$p.value))
cor$p=(nrow(cor):1)
p1=ggplot(data = cor)+
   geom_point(aes(x=cor,
                 y=p),
             alpha=0.5,
             size=1,
             color="gray")+
  #geom_text_repel(aes(x=cor,y=p,label=label))+
   theme_classic()+
   geom_vline(xintercept = 0,lty=5)+
   xlab("Correlation with granzyme expression")+
   ylab("Ranking")+
   annotate("rect",
            xmin=min(cor[1:top,"cor"]),
            xmax=Inf,
            ymin = 17000, ymax = Inf,
            alpha=0.5,
            fill="orange")


data = cor[1:top,]
data$label=factor(data$label,levels = rev(data$label))
p2=ggplot(data = data)+
  geom_segment(aes(xend=cor,
                   x=0,
                   y=label,
                   yend=label),
               color="grey") +
  geom_point(aes(x=cor,
                 y=label,
                 color=p.value),
             size=4)+
  scale_color_gradient(high = "#FFE4B5",low = "orange")+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("R") +
  ylab(" ")+
  labs(title = paste("Top Gene: ",top))

p1+p2+plot_layout(widths = c(2,1))

ggsave(filename = "Fig/Fig5e3.pdf",height = 6,width = 10)


