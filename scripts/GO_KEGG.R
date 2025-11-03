library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(enrichplot)
library(data.table)

rm(list=ls())
load("scRNA_endo.rds")

################################################################################
meta=scRNA@meta.data

Idents(scRNA)<-"CellFromTumor"

deg <-FindMarkers(scRNA,
                  ident.1 ="1",
                  ident.2 ="0")

#################################################################################################

sig_deg<- subset(deg, p_val_adj<0.01 & abs(avg_log2FC)>1)

ego_ALL <- enrichGO(gene          = row.names(sig_deg),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

ego_all <- data.frame(ego_ALL)

ego_BP <- enrichGO(gene          = row.names(sig_deg),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

dotplot(ego_ALL,split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scales = "free")

## KEGG富集分析
#基因ID转换
genelist <- bitr(row.names(sig_deg), 
                 fromType="SYMBOL",
                 toType="ENTREZID", 
                 OrgDb='org.Hs.eg.db')

genelist <- pull(genelist, ENTREZID)

ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')

dotplot(ekegg)

kegg <-data.frame(ekegg)


#  GSEA
###############################################################################

df <- bitr(unique(rownames(deg)), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG$gene =rownames(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='gene')

data_all_sort <- DEG %>% 
  arrange(desc(avg_log2FC))

geneList = data_all_sort$avg_log2FC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

gsea= gseKEGG(geneList = geneList)
dotplot(gsea)



## GSEA
###############################################################################
## fgsea中输入的关键基因信息
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- rownames(alldiff)
id
## fgsea中输入的关键通路信息
gmtfile <- "h.all.v7.4.symbols.gmt"

hallmark <- read.gmt(gmtfile)

hallmark$term <- gsub('HALLMARK_','',hallmark$term)

hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

## Perform the fgsea analysis
fgseaRes <- fgsea(pathways = hallmark.list, 
                  stats = id,
                  minSize=1,
                  maxSize=10000,
                  nperm=100)
sig <- fgseaRes
sig <- sig[order(sig$NES,decreasing = T),]
up=sig[1:10]$pathway
sig <- sig[order(sig$NES,decreasing = F),]
down=sig[1:10]$pathway
## 最后一步 开始绘图
dev.off()
plotGseaTable(hallmark.list[c(up,down)],id,fgseaRes,gseaParam = 0.5)

#################################################################################

alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]

id <- alldiff$avg_log2FC
names(id) <- rownames(alldiff)
id
gmtfile <- "h.all.v7.4.symbols.gmt"

hallmark <- read.gmt(gmtfile)

y <- GSEA(id,TERM2GENE = hallmark)

dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)

yd <- data.frame(y)

gseaplot2(y,"HALLMARK_HYPOXIA",color = "red",pvalue_table = T)
