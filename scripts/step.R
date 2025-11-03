library(data.table)
library(tidyverse)

LUAD <-fread("TCGA-LUAD.htseq_counts.tsv.gz",header = T,data.table = F)
LUSC <-fread("TCGA-LUSC.htseq_counts.tsv.gz",header = T,data.table = F)
plat <- data.table::fread(file = "gencode.v22.annotation.gene.probeMap",data.table = F)
str(LUAD)
str(LUSC)
colnames(plat)[1] <- c("Ensembl_ID")
plat <- plat[,1:2]
colnames(plat)

LUAD <- LUAD %>% 
  inner_join(plat,by="Ensembl_ID") %>% 
  select(Ensembl_ID,gene,everything()) %>% 
  select(-Ensembl_ID) %>%
  mutate(rowMean =rowMeans(.[grep("TCGA", names(.))])) %>% #求出平均数
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(gene,.keep_all = T) %>% # gene留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames("gene")


LUSC <- LUSC %>% 
  inner_join(plat,by="Ensembl_ID") %>% 
  select(Ensembl_ID,gene,everything()) %>% 
  select(-Ensembl_ID) %>%
  mutate(rowMean =rowMeans(.[grep("TCGA", names(.))])) %>% #求出平均数
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(gene,.keep_all = T) %>% # gene留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames("gene")

LUAD.id <-substring(colnames(LUAD),14,15)
table(LUAD.id)

LUSC.id <-substring(colnames(LUSC),14,15)
table(LUSC.id)

rev_log2 <- function(data){
  2^data - 1
}

LUAD <- rev_log2(LUAD)

LUSC <- rev_log2(LUSC)

library(IOBR)

LUAD_tpm <-count2tpm(LUAD,idType ="SYMBOL")
dim(LUAD_tpm)
LUAD_tpm <-log2(LUAD_tpm+1)

LUSC_tpm <-count2tpm(LUSC,idType ="SYMBOL")
dim(LUSC_tpm)
LUSC_tpm <-log2(LUSC_tpm+1)

lung1 <-LUAD_tpm[,which(LUAD.id=="11")]
dim(lung1)
lung2 <-LUSC_tpm[,which(LUSC.id=="11")]
dim(lung2)

LUAD_tumor <- LUAD_tpm[,which(LUAD.id!="11")]
LUSC_tumor <- LUSC_tpm[,which(LUSC.id!="11")]
## lung (n= 108), LUSC (n= 501) or LUAD (n= 513)

lung=cbind(lung1,lung2)

save(lung,lung1,lung2,LUAD_tumor,LUSC_tumor,file = "TCGA.Rds")
save(LUAD_tpm,LUSC_tpm,file = "TCGA_all.Rds")


























