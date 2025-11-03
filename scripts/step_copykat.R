library(devtools)
install_github("navinlabcode/copykat")
library(copykat)
library(Seurat)
load("scRNA_endo.Rds")

DimPlot(scRNA)

expr <- as.matrix(scRNA@assays$RNA@counts)

copykat.test <- copykat(rawmat=expr, 
                        id.type="S", 
                        cell.line="no", 
                        ngene.chr=3, 
                        LOW.DR = 0.01,
                        win.size=25, 
                        KS.cut=0.01, 
                        sam.name="st", 
                        distance="euclidean", 
                        n.cores=1)
save(copykat.test,file = "copy.Rds")

pre <-copykat.test$prediction
pre <-as.data.frame(pre)
pre <-pre[rownames(scRNA@meta.data),]

scRNA$copykat <-pre$copykat.pred

DimPlot(scRNA,group.by = "copykat")
