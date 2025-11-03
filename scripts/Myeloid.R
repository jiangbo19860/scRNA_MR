rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)

load("Myeloid.Cellview.Rds")
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
scale.genes <-  VariableFeatures(scRNA)    #内存不够可用此命令代替上一行
scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 

DimPlot(scRNA, reduction = "pca", group.by="orig.ident")

ElbowPlot(scRNA, ndims=20, reduction="pca") 

pc.num=1:15

scRNA <- FindNeighbors(scRNA, dims = pc.num) 

scRNA <- FindClusters(scRNA, resolution = 0.4)

metadata <- scRNA@meta.data

scRNA <- RunTSNE(scRNA, dims = pc.num)
scRNA <- RunUMAP(scRNA, dims = pc.num)

save(scRNA,file = "scRNA_Myeloid.Rds")

load("scRNA_Myeloid.Rds")

DimPlot(scRNA,group.by = "dbCluster",reduction = "tsne")
DimPlot(scRNA,reduction = "tsne")

DimPlot(scRNA,reduction = "tsne",label = T)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "none")+
  labs(title = "Myeloid")

ggsave(filename = "Fig/Fig4e.pdf",height = 8,width = 8)

####################################################################################

DimPlot(scRNA, 
        reduction = "tsne",
        group.by ="CellFromTumor")+
  scale_color_manual(values = c("#72C667","#2D5474"),
                     labels=c("Non-malignant","Tumor"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")

ggsave(filename = "Fig/Fig4e1.pdf",height = 8,width = 8)


####################################################################################
tmp=scRNA@meta.data
tmp$group=ifelse(tmp$CellFromTumor=="1","Tumor","Non-malignant")
tmp$seurat_clusters =factor(tmp$seurat_clusters,levels = 11:0)

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

ggsave("Fig/Fig4_Myeloid_bar.pdf",height = 4,width = 10)


################################################################################

p=FeaturePlot(scRNA, 
              reduction = "tsne",
              features = c("FCER1A",
                           "CD163",
                           "IFITM3",
                           "CLEC9A",
                           "S100A12",
                           "RGCC"),
              ncol = 3,
              cols = c("gray","red"))
p & theme(panel.border = element_rect(fill=NA,color="black", 
                                      size=1, linetype="solid"),
          legend.position = "none")
p
ggsave("Fig/Fig4f.pdf",height = 6,width = 9)

######################################################################################

DimPlot(scRNA,reduction = "tsne",label = T)+
  FeaturePlot(scRNA, 
              reduction = "tsne",
              features = c("CD163"),
              cols = c("gray","red"),
              label = T)

library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

gs <-read.gmt("h.all.v7.4.symbols.gmt")

gs$term <- gsub('HALLMARK_','',gs$term)

gs.list <- gs %>% split(.$term) %>% lapply( "[[", 2)

scRNAsub <-subset(scRNA,seurat_clusters %in% c(0,1,3,5:8,10,11))

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

head(x)

write.csv(x, "Myeloid_gsva_limma.csv", quote = F)

df <- data.frame(ID = rownames(x), score = x$t)

#按照score的值分组
cutoff <- 4
df$group <- cut(df$score, 
                breaks = c(-Inf, -cutoff, cutoff, Inf),
                labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, 
                   score, 
                   fill = group)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 
                               'snow3', 
                               'dodgerblue4'), 
                    guide = "none") + 
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + 
  geom_text(data = subset(df, score < 0),
            aes(x=ID, 
                y= 0.1, 
                label= ID, 
                color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = 0 ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, 
                y= -0.1, 
                label=ID, 
                color = group),
            size = 3, 
            hjust = 1) +  
  scale_colour_manual(values = c("black","snow3","black"), 
                      guide = "none") +
  xlab("") +
  ylab("t value of GSVA score, tumor \n versus non-malignant")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) #去除y轴

ggsave(filename = "Fig/Fig4h.pdf",height = 12,width = 10)

#################################################################################
