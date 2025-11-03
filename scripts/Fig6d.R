library(ggpubr)
library(Seurat)
library(ggplot2)
library(survival)
library(ggthemes)
library(tidyverse)

load("TCGA/TCGA.Rds")
##rt作为输入 第一二列分别为 status 和 survival 
# 名称一定是这两个 后面列是细胞亚群
load("scRNA.Rds")
Idents(scRNA)<-"ClusterName"
DimPlot(scRNA)
dput(unique(scRNA$ClusterName))

cells <-c("CD4+ T cells ", "CD8+ T cells ", "mast cells ", "Langerhans cells ", 
          "lower quality endothelial cell", "macrophages", "cancer cells pt 4", 
          "regulatory T cells ", "plasma B cells ", "tumour endothelial cell")

scRNA <-subset(scRNA,ClusterName %in% cells)

sc.marker <-FindAllMarkers(scRNA,
                           only.pos = T,
                           min.pct = 0.25)

phe_LUAD = data.table::fread("TCGA/TCGA-LUAD.survival.tsv",header = T,data.table = F)
rownames(phe_LUAD) =phe_LUAD$sample
phe_LUAD =phe_LUAD[colnames(LUAD_tumor),]
rownames(phe_LUAD)<-colnames(LUAD_tumor)

phe_LUSC = data.table::fread("TCGA/TCGA-LUSC.survival.tsv",header = T,data.table = F)
rownames(phe_LUSC) =phe_LUSC$sample
phe_LUSC =phe_LUSC[colnames(LUSC_tumor),]
rownames(phe_LUSC)<-colnames(LUSC_tumor)

data <-data.frame(row.names = colnames(LUAD_tumor))

for (x in unique(sc.marker$cluster)) {
  marker_tmp <-subset(sc.marker,cluster==x)
  gene <-intersect(marker_tmp$gene,rownames(LUAD_tumor))
  data_tmp <-as.data.frame(LUAD_tumor[gene,])
  
  for (i in 1:ncol(data_tmp)){
    data_tmp["avr",i]<-mean(as.numeric(data_tmp[1:length(gene),i]))
  }
  
  data[,x]<-as.numeric(data_tmp["avr",])
}

rt<-data.frame(row.names = colnames(LUAD_tumor),
               survival=phe_LUAD$OS.time,
               status=phe_LUAD$OS)
rt<-cbind(rt,data)

outTab=data.frame()

for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(survival, status) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                            z=coxSummary$coefficients[,"z"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

write.csv(outTab,"outTab.csv")
outTab=read.csv("outTab.csv",header = T)
###这里为什么要写出再读入？因为这样避免结果的数值被认为是字符串
outTab$group=ifelse(outTab$pvalue>0.05,"a",
                    ifelse(outTab$z>0,"b","c"))

ggplot(outTab, aes(x=gene, y=z,color=group)) + 
  geom_point(stat='identity',size=3)  +
  scale_color_manual(values=c("gray","blue","lightblue"))+
  geom_segment(aes(y = 0, 
                   x = gene, 
                   yend = z, 
                   xend = gene), 
               color = "grey",
               linetype="dashed") +
  scale_y_continuous(limits = c(-5,5))+   ##x轴长度 如果下面有warning一般是x轴小于z值
  theme_base() +
  theme(legend.position ="")+
  coord_flip() +
  geom_hline(yintercept = c(2,-2),linetype="dashed",colour="grey",size=1)+ #虚线
  xlab("")+ylab("Z score")+
  geom_hline(yintercept = c(3,-3),linetype="solid")+  #实线
  labs(title="LUAD")

####################################################################################

library(survminer)
data_sur=rt
for (j in 3:ncol(data_sur)) {
  data_sur[,j]<-ifelse(as.numeric(data_sur[,j])>median(data_sur[,j]),
                       "high",
                       "low")
}

plot<-list()
for (j in 1:10) {
  data_sur_plot <-data_sur[,c(1,2,j+2)]
  colnames(data_sur_plot)<-c("survival","status","group")
  km <-ggsurvplot(survfit(Surv(survival, status) ~ group,
                          data_sur_plot),
                  pval = T)
  plot[j]<-km[1]
}
p.all <- cowplot::plot_grid(plotlist = plot,ncol = 2)
save_plot("Fig/Fig6d.pdf", p.all,base_height = 20,base_width = 10)




