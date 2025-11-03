rm(list = ls())
load("All_cell_subcluster.Rds")
load("region.Rds")
meta=subset(meta2,region!="middle" & cell !="Cancer")

table(meta$PatientNumber,meta$region)

meta$Cluster=All.meta[rownames(meta),"dbCluster"]

meta$cell=All.meta[rownames(meta),"cell"]

meta$sub=paste0(meta$cell,"_",meta$Cluster)

meta$region2=paste0(meta$region,"_",meta$PatientNumber)

data_sub=as.data.frame(table(meta$sub,meta$region2))

data=as.data.frame(table(meta$PatientNumber,meta$region))

# packages
library(reshape2)
library(tidyr)

data_plot=spread(data_sub,Var1,Freq)
rownames(data_plot)=data_plot$Var2
data_plot=data_plot[,-1]
data_plot$total=data$Freq

for (i in 1:ncol(data_plot)) {
  data_plot[,i]=as.numeric(data_plot[,i])/(data_plot$total)
}

data_plot=data_plot[,-ncol(data_plot)]
data_plot$group=as.numeric(c(rep(0,5),rep(1,5)))

dput(unique(meta$sub))

cell_name=c("B_cell_1", "EC_2", "B_cell_4", "B_cell_5", "B_cell_3", "T_cell_1", 
            "Fibroblast_5", "B_cell_2", "Myeloid_3", "T_cell_3", "Fibroblast_3", 
            "T_cell_9", "T_cell_5", "T_cell_2", "B_cell_7", "B_cell_6", "T_cell_6", 
            "EC_4", "Myeloid_2", "Myeloid_5", "Myeloid_11", "Myeloid_12", 
            "T_cell_8", "Alveolar_1", "T_cell_7", "Fibroblast_2", "B_cell_8", 
            "EC_6", "Alveolar_2", "Myeloid_9", "Myeloid_7", "Myeloid_1", 
            "EC_3", "Epithelial_2", "Alveolar_4", "Fibroblast_6", "Alveolar_7", 
            "T_cell_4", "Alveolar_3", "Fibroblast_4", "Myeloid_10", "Myeloid_4", 
            "Fibroblast_1", "Myeloid_6", "Epithelial_1", "Myeloid_8", "EC_5", 
            "Alveolar_5", "EC_1", "Fibroblast_7")

data_plot2=data.frame(cell=NULL,
                      t=NULL,
                      p=NULL)

for (i in cell_name) {
  data_tmp=data_plot[,c(i,"group")]
  x=lm(group~.,data = data_tmp)
  t=summary(x)$coefficients[,3] 
  p=summary(x)$coefficients[,4]
  data_tmp2=data.frame(cell=i,
                       t=t[2],
                       p=p[2])
  data_plot2=rbind(data_plot2,data_tmp2)
}


library(ggpubr)

data_plot2$pvalue=ifelse(data_plot2$p>0.1,"c",ifelse(data_plot2$p>0.05,"b","a"))
data_plot2$major=data.frame(sapply(data_plot2$cell, function(x)unlist(strsplit(x,'_'))[1]),stringsAsFactors = F)[,1]


library(cowplot)

ggplot(data_plot2, aes(x=cell, y=t, fill=pvalue,shape=pvalue))+ 
  geom_segment(aes(y = 0, 
                   x = cell, 
                   yend = t, 
                   xend = cell), 
               color = "#36648B")+
  geom_point(stat='identity', 
             size =3, shape = 21, 
             stroke = 1,color="#36648B")  +
  scale_fill_manual(values = c("#000080","#4F94CD","#FFFFFF"),
                    name="P value",
                    breaks=c("a","b","c"),
                    labels=c("p <0.05","p<0.1","P>0.1"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  scale_x_discrete(position = "top")+
  scale_y_reverse()+
  geom_hline(yintercept =0,lty=5)+
  geom_hline(yintercept =1.7,lty=5,color="gray")+
  geom_hline(yintercept =2.2,lty=1,color="gray")

ggplot(data_plot2, aes(x=cell, y=t, fill=pvalue,shape=pvalue))+ 
  geom_segment(aes(y = 0, 
                   x = cell, 
                   yend = t, 
                   xend = cell), 
               color = "#36648B")+
  geom_point(stat='identity', 
             size =3, shape = 21, 
             stroke = 1,color="#36648B")  +
  scale_fill_manual(values = c("#000080","#4F94CD","#FFFFFF"),
                    name="P value",
                    breaks=c("a","b","c"),
                    labels=c("p <0.05","p<0.1","P>0.1"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  scale_x_discrete(position = "top")+
  scale_y_reverse()+
  geom_hline(yintercept =0,lty=5)+
  geom_hline(yintercept =1.7,lty=5,color="gray")+
  geom_hline(yintercept =2.2,lty=1,color="gray")+
  annotate("rect",
           xmin="B_cell_1",
           xmax="B_cell_8",
           ymin = -Inf, 
           ymax = Inf,
           alpha=0.5,
           fill="gray")+
  annotate("rect",
           xmin="Epithelial_1",
           xmax="Epithelial_2",
           ymin = -Inf, 
           ymax = Inf,
           alpha=0.5,
           fill="gray")+
  annotate("rect",
           xmin="Myeloid_1",
           xmax="Myeloid_2",
           ymin = -Inf, 
           ymax = Inf,
           alpha=0.5,
           fill="gray")+
  annotate("rect",
           xmin="Myeloid_1",
           xmax="Myeloid_9",
           ymin = -Inf, 
           ymax = Inf,
           alpha=0.5,
           fill="gray")

ggsave(filename = "Fig/Fig6b.pdf",height = 6,width = 10)










