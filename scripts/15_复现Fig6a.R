rm(list = ls())

load("scRNA.Rds")
load("Allsamples.Cellview.Rds")

metadata<-scRNA@meta.data

meta2=subset(metadata, CellFromTumor=="1")

table(meta2$PatientNumber,meta2$Sample)

# c(1,10,14,19,21) "middle"
# c(4,12,13,18,22) "edge"
# c(5,8,15,20,23) "core"

meta2$region=ifelse(meta2$Sample %in% c(1,10,14,19,21),"middle",
                    ifelse(meta2$Sample %in% c(4,12,13,18,22),"edge","core"))

table(meta2$PatientNumber,meta2$region)

# core edge middle
# 1  114 1651    806
# 2  424  648    400
# 3 2149 3554   3266
# 4 5270 3283   5198
# 5 4913 4664   2983

meta2$cell=meta2$ClusterName

meta2[grep("cancer cells",meta2$ClusterName),]$cell="Cancer"

meta3=subset(meta2,cell=="Cancer")

table(meta3$PatientNumber,meta3$region)

## 计算肿瘤细胞比例，首先要计算细胞总量

pop1=as.data.frame(table(meta2$PatientNumber,meta2$region))

## 提取肿瘤细胞数量
pop=as.data.frame(table(meta3$PatientNumber,meta3$region))

pop$all=pop1$Freq
## 计算比例
pop$p=(pop$Freq)/(pop$all)

pop$Var2=factor(pop$Var2,levels = c("core","middle","edge"))

top=ggbarplot(data = pop,
              x="Var1",
              y="p",
              color = "Var2",
              fill = "Var2",
              position = position_dodge(0.8),
              palette = c("#284B66","#4F94CC","#78C1FF"),
              ylab = "Cancer",
              xlab = " ")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

cancer=scRNA[,rownames(meta3)]
cancer=GetAssayData(cancer)
cancer=as.data.frame(cancer)  
hypoxia_marker <-read.csv("hypoxia.csv",header = T)
gene <-intersect(hypoxia_marker$Symbol,rownames(cancer))

cancer <-cancer[gene,]

for (i in 1:ncol(cancer)) {
  cancer["hypoxia",i]<-mean(as.numeric(cancer[1:46,i]))
}

data=data.frame(hypoxia=as.numeric(cancer["hypoxia",]),
                row.names = colnames(cancer))
data=cbind(data,meta3[,c("region","PatientNumber")])
data$region=factor(data$region,levels = c("core","middle","edge"))

bottom=ggbarplot(data = data,
                 x="PatientNumber",
                 y="hypoxia",
                 color = "region",
                 fill = "region",
                 position = position_dodge(0.8),
                 palette = c("#284B66","#4F94CC","#78C1FF"))+
  theme(legend.position  = "right",
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

top+bottom+plot_layout(ncol = 1,guides = "collect")

save(meta2,meta3,file = "region.Rds")


