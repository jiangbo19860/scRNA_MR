rm(list = ls())
library(Seurat)
library(ggplot2)

load("TCGA/TCGA.Rds")
load("scRNA_endo.Rds")

FeaturePlot(scRNA,"MT2A",label = T,reduction = "tsne")  #normal
FeaturePlot(scRNA,"HSPG2",label = T,reduction = "tsne") #tumor

Idents(scRNA)<-"ClusterName"
DimPlot(scRNA,label = T)

marker=FindMarkers(scRNA,
                   ident.1 = "normal endothelial cell",
                   ident.2 = "tumour endothelial cell")

Normal_EC <-rownames(marker[order(marker$avg_log2FC,decreasing = T),])[1:10]

Tumor_EC  <-rownames(marker[order(marker$avg_log2FC,decreasing = F),])[1:10]


data_lung <-as.data.frame(lung[Normal_EC,])

for (i in 1:ncol(data_lung)){
  data_lung["avr",i]<-mean(as.numeric(data_lung[1:10,i]))
}

data_LUAD <-as.data.frame(LUAD_tumor[Normal_EC,])

for (i in 1:ncol(data_LUAD)){
  data_LUAD["avr",i]<-mean(as.numeric(data_LUAD[1:10,i]))
}

data_LUSC <-as.data.frame(LUSC_tumor[Normal_EC,])

for (i in 1:ncol(data_LUSC)){
  data_LUSC["avr",i]<-mean(as.numeric(data_LUSC[1:10,i]))
}

data_plot <- data.frame(Avr=c(as.numeric(data_lung["avr",]),
                              as.numeric(data_LUSC["avr",]),
                              as.numeric(data_LUAD["avr",])),
                       type=c(rep("Lung",108),
                              rep("LUSC",501),
                              rep("LUAD",526)))

ggpubr::ggboxplot(data = data_plot,
                  x="type",
                  y="Avr",
              fill ="type",
              palette = c("#7CCD7C",
                          "#30597C",
                          "#5A91BF"),
              xlab = "",
              ylab = "Average gene expression")+
  theme_bw()+
  theme(legend.position = "none")

###########################################################################################

data_lung <-as.data.frame(lung[Tumor_EC,])

for (i in 1:ncol(data_lung)){
  data_lung["avr",i]<-mean(as.numeric(data_lung[1:10,i]))
}

data_LUAD <-as.data.frame(LUAD_tumor[Tumor_EC,])

for (i in 1:ncol(data_LUAD)){
  data_LUAD["avr",i]<-mean(as.numeric(data_LUAD[1:10,i]))
}

data_LUSC <-as.data.frame(LUSC_tumor[Tumor_EC,])

for (i in 1:ncol(data_LUSC)){
  data_LUSC["avr",i]<-mean(as.numeric(data_LUSC[1:10,i]))
}

data_plot <- data.frame(Avr=c(as.numeric(data_lung["avr",]),
                              as.numeric(data_LUSC["avr",]),
                              as.numeric(data_LUAD["avr",])),
                        type=c(rep("Lung",108),
                               rep("LUSC",501),
                               rep("LUAD",526)))

ggpubr::ggboxplot(data = data_plot,
                  x="type",
                  y="Avr",
                  fill ="type",
                  palette = c("#7CCD7C",
                              "#30597C",
                              "#5A91BF"),
                  xlab = "",
                  ylab = "Average gene expression")+
  theme_bw()+
  theme(legend.position = "none")
######################################################################################

## IGFBP3, HSPG2, INSR, PLVAP, ENPP2

#######################################################################################
Tumor_EC <-c("IGFBP3", "HSPG2", "INSR", "PLVAP", "ENPP2")

data_lung <-as.data.frame(lung[Tumor_EC,])

for (i in 1:ncol(data_lung)){
  data_lung["avr",i]<-mean(as.numeric(data_lung[1:5,i]))
}

data_LUAD <-as.data.frame(LUAD_tumor[Tumor_EC,])

for (i in 1:ncol(data_LUAD)){
  data_LUAD["avr",i]<-mean(as.numeric(data_LUAD[1:5,i]))
}

data_LUSC <-as.data.frame(LUSC_tumor[Tumor_EC,])

for (i in 1:ncol(data_LUSC)){
  data_LUSC["avr",i]<-mean(as.numeric(data_LUSC[1:5,i]))
}

data_plot <- data.frame(Avr=c(as.numeric(data_lung["avr",]),
                              as.numeric(data_LUSC["avr",]),
                              as.numeric(data_LUAD["avr",])),
                        type=c(rep("Lung",108),
                               rep("LUSC",501),
                               rep("LUAD",526)))

ggpubr::ggboxplot(data = data_plot,
                  x="type",
                  y="Avr",
                  fill ="type",
                  palette = c("#7CCD7C",
                              "#30597C",
                              "#5A91BF"),
                  xlab = "",
                  ylab = "Average gene expression")+
  theme_bw()+
  theme(legend.position = "none")

########################################################################################

#

#
#   Fibroblast
#

#

#######################################################################################

rm(list = ls())
library(Seurat)
library(ggplot2)
library(tidyverse)

load("TCGA/TCGA.Rds")
load("scRNA_fibro.Rds")

DimPlot(scRNA,label = T)

marker=FindAllMarkers(scRNA)

Cluster <-marker %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)

Cluster0 <-Cluster[which(Cluster$cluster=="0"),]$gene

data_lung <-as.data.frame(lung[Cluster0,])

for (i in 1:ncol(data_lung)){
  data_lung["avr",i]<-mean(as.numeric(data_lung[1:10,i]))
}

data_LUAD <-as.data.frame(LUAD_tumor[Cluster0,])

for (i in 1:ncol(data_LUAD)){
  data_LUAD["avr",i]<-mean(as.numeric(data_LUAD[1:10,i]))
}

data_LUSC <-as.data.frame(LUSC_tumor[Cluster0,])

for (i in 1:ncol(data_LUSC)){
  data_LUSC["avr",i]<-mean(as.numeric(data_LUSC[1:10,i]))
}

data_plot <- data.frame(Avr=c(as.numeric(data_lung["avr",]),
                              as.numeric(data_LUSC["avr",]),
                              as.numeric(data_LUAD["avr",])),
                        type=c(rep("Lung",108),
                               rep("LUSC",501),
                               rep("LUAD",526)))

ggpubr::ggboxplot(data = data_plot,
                  x="type",
                  y="Avr",
                  fill ="type",
                  palette = c("#7CCD7C",
                              "#30597C",
                              "#5A91BF"),
                  xlab = "",
                  ylab = "Average gene expression")+
  theme_bw()+
  theme(legend.position = "none")

###########################################################################################








