# devtools::install_github("sqjin/CellChat")
rm(list = ls())
load("scRNA_endo.rds")
library(CellChat)
library(Seurat)
options(stringsAsFactors = FALSE)

data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data") # normalized data matrix
Idents(scRNA)<-"ClusterName"
labels <- Idents(scRNA)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "group")

cellchat <- addMeta(cellchat, meta = meta)

cellchat <- setIdent(cellchat, ident.use = "group") # set "labels" as default cell identity

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multiprocess", workers = 4) # do scRNArallel

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
dev.off()
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= T, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize,
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
dev.off()
library(tidyverse)
df.net <- subsetCommunication(cellchat,thresh = 0.05)
# "normal endothelial cell"       
# "tumour endothelial cell"
df<-df.net%>%filter(source=="tumour endothelial cell")
df<-df.net%>%filter(target=="tumour endothelial cell")

dev.off()
pathways.show <- c("MIF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, 
                    signaling =pathways.show,  
                    vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.


dev.off()
unique(labels)
netVisual_bubble(cellchat, 
                 sources.use = 15, 
                 targets.use = 6, 
                 remove.isolate = FALSE)

netVisual_bubble(cellchat, 
                 sources.use = 2, 
                 targets.use = c(1:7), 
                 remove.isolate = FALSE)

netVisual_chord_gene(cellchat, 
                     sources.use = c(1:7), 
                     targets.use =2, 
                     lab.cex =0.3,
                     legend.pos.y = 20)


StackedVlnPlot(scRNA,features = c("EPCAM","PTPRC"))

load("position.Rds")
load("stRNA.Rds")
embed_umap2=data.frame(UMAP_1=position_sub_sub$x,
                       UMAP_2=position_sub_sub$y,
                       row.names = rownames(position_sub_sub))

stRNA@reductions$umap@cell.embeddings=as.matrix(embed_umap2)
DimPlot(stRNA)
FeaturePlot(stRNA,"CD99")

saveRDS(cellchat, file = "cellchat.rds")

cellchat=readRDS("cellchat.rds")

