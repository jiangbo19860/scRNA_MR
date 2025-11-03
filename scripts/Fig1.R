load("scRNA.Rds")
# Fig1b
dir.create("Fig")
metadata<-scRNA@meta.data

DimPlot(scRNA, 
        reduction = "tsne",
        group.by ="CellFromTumor")+
  scale_color_manual(values = c("#72C667","#2D5474"),
                     labels=c("Non-malignant","Tumor"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")

ggsave(filename = "Fig/Fig1b.pdf",height = 8,width = 8)

#####################################################################################
DimPlot(scRNA, 
        reduction = "tsne",
        group.by ="PatientNumber")+
  scale_color_manual(values = c("#268A24",
                                "#8EEE8B",
                                "#FB6346",
                                "#FBD51A",
                                "#28507D"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Patient")

ggsave(filename = "Fig/Fig1b2.pdf",height = 8,width = 8)

#####################################################################################

DimPlot(scRNA, 
        reduction = "tsne",
        group.by ="dbCluster",label = T)+
  scale_color_manual(values = c("#3A6135",
                                "#E7A649",
                                "#4C6679",
                                "#814ECC",
                                "#6466A5",
                                "#A14299",
                                "#F58C70",
                                "#B4434E"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "none")+
  labs(title = "Cell type")


metadata$cell=metadata$ClusterName

metadata[grep("T cells|natural killer",metadata$ClusterName),]$cell="T cells"
metadata[grep("B cells|granulocytes",metadata$ClusterName),]$cell="B cells"
metadata[grep("fibroblasts",metadata$ClusterName),]$cell="Fibroblasts"
metadata[grep("endothelial cell",metadata$ClusterName),]$cell="Endothelial"
metadata[grep("epithelial cell|EC|basal cells",metadata$ClusterName),]$cell="Epithelial"
metadata[grep("cancer cells",metadata$ClusterName),]$cell="Cancer"
metadata[grep("alveolar",metadata$ClusterName),]$cell="Alveolar"
metadata[grep("macrophages|Langerhans|mast cells|dendritic",metadata$ClusterName),]$cell="Myleoid"
metadata[grep("erythroblasts|secretory club cells",metadata$ClusterName),]$cell=NA

table(metadata$cell)
scRNA$cell=metadata$cell

DimPlot(scRNA, reduction = "tsne",
        group.by ="cell",label = T)+
  scale_color_manual(values = c("#3A6135",
                                "#E7A649",
                                "#4C6679",
                                "#814ECC",
                                "#6466A5",
                                "#A14299",
                                "#F58C70",
                                "#B4434E"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "none")+
  labs(title = "Cell type")

ggsave(filename = "Fig/Fig1b3.pdf",height = 8,width = 8)

####################################################################################

FeaturePlot(scRNA, 
            reduction = "tsne",
            features = "nUMI")+
  scale_color_gradient(high = "blue",low = "gray")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "none")+
  labs(title = "Transcript count")

ggsave(filename = "Fig/Fig1b4.pdf",height = 8,width = 8)


pdf("Fig/Fig1c.pdf",height = 6,width = 12)
p=FeaturePlot(scRNA, reduction = "tsne",
              features = c("CLDN18",
                           "CLDN5",
                           "CAPS",
                           "COL1A1",
                           "CD79A",
                           "LYZ",
                           "CD3D",
                           "EPCAM"),
              ncol = 4,
              cols = c("gray","red"))
p & theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          legend.position = "none")
dev.off()

######################################################################################


