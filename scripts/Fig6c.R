rm(list = ls())
library(Seurat)
library(tidyverse)

load("scRNA.Rds")
metadata <-scRNA@meta.data
metadata$cell=metadata$ClusterName

metadata[grep("T cells|natural killer",metadata$ClusterName),]$cell="T cells"
metadata[grep("B cells|granulocytes",metadata$ClusterName),]$cell="B cells"
metadata[grep("fibroblasts",metadata$ClusterName),]$cell="Fibroblasts"
metadata[grep("endothelial cell",metadata$ClusterName),]$cell="Endothelial"
metadata[grep("epithelial cell|EC|basal cells",metadata$ClusterName),]$cell="Epithelial"
metadata[grep("cancer cells",metadata$ClusterName),]$cell="Cancer"
metadata[grep("alveolar",metadata$ClusterName),]$cell="Alveolar"
metadata[grep("macrophages|Langerhans|mast cells|dendritic",metadata$ClusterName),]$cell="Myleoid"
metadata[grep("erythroblasts|secretory club cells",metadata$ClusterName),]$cell="basel"

metadata -> scRNA@meta.data
Idents(scRNA)<-"cell"
DimPlot(scRNA)

# sc.marker <-FindAllMarkers(scRNA,
#                            only.pos = T,
#                            logfc.threshold = 0.5,
#                            min.pct = 0.25)

load("scRNA_fibro.Rds")

sc.marker <-FindAllMarkers(scRNA,only.pos = T,logfc.threshold = 0.5,min.pct = 0.25)
sc.marker <-sc.marker %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)

load("TCGA/TCGA.Rds")

data_plot_all <- data.frame(LUAD=NULL,
                            LUSC=NULL,
                            fc=NULL,
                            p=NULL,
                         cluster=NULL)

for (x in 0:5) {
  sc.marker.tmp <- subset(sc.marker, cluster == x)
  marker <- intersect(rownames(lung), sc.marker.tmp$gene)
  
  data_LUAD <- as.data.frame(LUAD_tumor[marker, ])
  for (i in 1:ncol(data_LUAD)) {
    data_LUAD["avr", i] <- mean(as.numeric(data_LUAD[1:length(marker), i]))
  }
  
  data_LUSC <- as.data.frame(LUSC_tumor[marker, ])
  for (i in 1:ncol(data_LUSC)) {
    data_LUSC["avr", i] <- mean(as.numeric(data_LUSC[1:length(marker), i]))
  }
  
  a = as.numeric(data_LUSC["avr", ])
  b = as.numeric(data_LUAD["avr", ])
  p  = wilcox.test(a, b, na.rm = T)$p.value
  fc = mean(a) / mean(b)
  data_plot <- data.frame(
    LUSC = mean(a),
    LUAD = mean(b),
    fc = fc,
    p = p,
    cluster = x
  )
  data_plot_all <- rbind(data_plot, data_plot_all)
}

data_plot_all$padj=p.adjust(data_plot_all$p,method = "fdr")

scale_lim=c((min(data_plot_all[,c("LUAD","LUSC")])-0.05),
            (max(data_plot_all[,c("LUAD","LUSC")])+0.05))

library(ggrepel)
data_plot_all$FC=ifelse(data_plot_all$fc>1,">1","<1")
data_plot_all$label=paste0("C",data_plot_all$cluster)

ggplot(data = data_plot_all)+
  geom_point(aes(x=LUSC,
                 y=LUAD,
                 color=FC),
             alpha=0.5)+
  scale_color_manual(values = c("black","red"))+
  geom_text_repel(aes(x=LUSC,
                      y=LUAD,
                      label=label))+
  theme_bw()+
  xlim(scale_lim)+ylim(scale_lim)+
  geom_abline(slope=1, intercept=0,lty=2)

