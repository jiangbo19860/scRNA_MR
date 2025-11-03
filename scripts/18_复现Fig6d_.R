# 复现Fig6d：LUAD和LUSC中单细胞亚群标志基因的Cox回归分析。整合了单细胞数据分析与临床生存分析

rm(list = ls()) # 清空工作空间
# 加载所需R包：
# ggpubr：增强绘图功能；Seurat：单细胞数据分析；ggplot2：可视化；
# survival：生存分析（Cox模型、KM曲线）；ggthemes：提供更多绘图主题；tidyverse：数据处理
pacman::p_load(
  ggpubr, Seurat, ggplot2, survival, ggthemes, tidyverse, here
)
here()

# 加载TCGA数据（包含LUAD_tumor、LUSC_tumor等肿瘤表达矩阵）
load("TCGA.Rds")

## 注释：后续生存分析的输入数据要求第一列为生存状态（status），第二列为生存时间（survival）
## 名称必须严格匹配，后续列将是各细胞亚群的标志基因平均表达量

# subset()筛选目标细胞亚群，FindAllMarkers()识别每个亚群marker genes(表达上调≥25%) --------
# 加载单细胞RNA-seq数据对象scRNA
load("scRNA.Rds")

# 设置scRNA的默认身份（Idents）为"ClusterName"（原始聚类名称）
Idents(scRNA) <- "ClusterName"

# 绘制降维图（UMAP/t-SNE），展示不同聚类（ClusterName）的分布
DimPlot(scRNA)

# 输出scRNA中所有独特的聚类名称（ClusterName），用于确认目标细胞亚群
dput(unique(scRNA$ClusterName))

# 定义需要分析的细胞亚群列表（从上述输出的ClusterName中筛选）
cells <- c("CD4+ T cells ", "CD8+ T cells ", "mast cells ", "Langerhans cells ", 
           "lower quality endothelial cell", "macrophages", "cancer cells pt 4", 
           "regulatory T cells ", "plasma B cells ", "tumour endothelial cell")

# 从scRNA中筛选出上述细胞亚群的数据（仅保留目标亚群进行后续分析）
scRNA <- subset(scRNA, ClusterName %in% cells)

# 计算每个细胞亚群的标志基因：
# only.pos=T：仅保留上调基因；min.pct=0.25：基因在至少25%的细胞中表达
sc.marker <- FindAllMarkers(scRNA,
                            only.pos = T,
                            min.pct = 0.25)


# 生存数据处理与表达量整合 ------------------------------------------------------------
# 读取LUAD的生存数据
# 格式：包含sample（样本ID）、OS（生存状态：0=存活，1=死亡）、OS.time（生存时间）
phe_LUAD <- data.table::fread("5_TCGA/TCGA-LUAD.survival.tsv", header = T, data.table = F)
# 确认列名和数据有效性
colnames(phe_LUAD) # [1] "sample"   "OS"       "_PATIENT" "OS.time" 
head(phe_LUAD) # 1 TCGA-NJ-A4YI-01A  1 TCGA-NJ-A4YI       4

# 将样本ID设置为行名（便于匹配表达数据）
rownames(phe_LUAD) <- phe_LUAD$sample
# 按LUAD_tumor的样本顺序筛选生存数据（确保一一对应）
phe_LUAD <- phe_LUAD[colnames(LUAD_tumor), ]
# 统一行名为LUAD_tumor的样本名
rownames(phe_LUAD) <- colnames(LUAD_tumor)

# 读取LUSC的生存数据
phe_LUSC <- data.table::fread("5_TCGA/TCGA-LUSC.survival.tsv", header = T, data.table = F)
# 确认列名和数据有效性
colnames(phe_LUSC)
head(phe_LUSC)

# 将样本ID设置为行名
rownames(phe_LUSC) <- phe_LUSC$sample
# 按LUSC_tumor的样本顺序筛选生存数据
phe_LUSC <- phe_LUSC[colnames(LUSC_tumor), ]
# 统一行名为LUSC_tumor的样本名
rownames(phe_LUSC) <- colnames(LUSC_tumor)

# 创建空数据框，存储LUAD样本中各细胞亚群标志基因的平均表达量
data <- data.frame(row.names = colnames(LUAD_tumor))

# 循环计算每个细胞亚群的标志基因在LUAD样本中的平均表达量
for (x in unique(sc.marker$cluster)) {  # x为亚群编号
  # 提取当前亚群的标志基因
  marker_tmp <- subset(sc.marker, cluster == x) #          p_val avg_log2FC pct.1 pct.2 p_val_adj                 cluster    gene SPARCL11     0   6.044733 0.945 0.018         0 tumour endothelial cell SPARCL1
  # 筛选出同时存在于LUAD_tumor和标志基因中的基因（确保匹配）
  gene <- intersect(marker_tmp$gene, rownames(LUAD_tumor))
  # 提取这些基因在LUAD_tumor中的表达数据
  data_tmp <- as.data.frame(LUAD_tumor[gene, , drop = FALSE])  # drop=FALSE避免单列时转为向量
  
  # 计算每个LUAD样本中当前亚群标志基因的平均表达量
  if (length(gene) > 0) {  # 仅当有有效基因时计算
    for (i in 1:ncol(data_tmp)) {
      data_tmp["avr", i] <- mean(as.numeric(data_tmp[1:length(gene), i]), na.rm = TRUE)
    }
    # 将平均表达量添加到data中
    data[, as.character(x)] <- as.numeric(data_tmp["avr", ])
  } else {
    # 若无匹配基因，填充NA并警告
    data[, as.character(x)] <- NA
    warning(paste("Cluster", x, "has no matching genes in LUAD_tumor"))
  }
}

head(marker_tmp)

# 构建生存分析输入数据框rt
rt <- data.frame(
  row.names = colnames(LUAD_tumor),
  survival = phe_LUAD$OS.time,  # 提取LUAD生存时间
  status = phe_LUAD$OS          # 提取LUAD生存状态（0/1）
)
rt <- cbind(rt, data)  # 合并生存信息与表达数据

# 数据有效性检查
cat("LUAD样本量：", nrow(rt), "\n")
cat("生存时间缺失情况：\n")
print(table(is.na(rt$survival)))
cat("生存状态缺失情况：\n")
print(table(is.na(rt$status)))
cat("样本名匹配数量：", length(intersect(colnames(LUAD_tumor), rownames(phe_LUAD))), "\n")
cat("前5行数据预览：\n")
print(head(rt))

# 移除全为NA的亚群列（避免Cox分析出错）
valid_cols <- colnames(rt)[!apply(rt, 2, function(col) all(is.na(col)))]
rt <- rt[, valid_cols, drop = FALSE]

# 创建空数据框存储Cox回归结果
outTab <- data.frame()

# 对每个细胞亚群进行Cox比例风险模型分析
if (ncol(rt) >= 3) {  # 确保有至少1个亚群
  for (i in colnames(rt[, 3:ncol(rt)])) {  # i为亚群名称
    # 提取当前亚群的非缺失数据
    valid_data <- rt[!is.na(rt[, i]) & !is.na(rt$survival) & !is.na(rt$status), ]
    
    if (nrow(valid_data) >= 10) {  # 样本量至少10才进行分析
      # 构建Cox模型
      cox <- coxph(Surv(survival, status) ~ valid_data[, i], data = valid_data)
      coxSummary <- summary(cox)
      # 存储结果
      outTab <- rbind(outTab, data.frame(
        gene = i,
        HR = coxSummary$coefficients[, "exp(coef)"],
        z = coxSummary$coefficients[, "z"],
        pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
      ))
    } else {
      warning(paste("Cluster", i, "has insufficient valid samples (<10), skipping Cox analysis"))
    }
  }
} else {
  warning("No valid cell clusters for Cox analysis")
}

# 处理Cox结果（避免数据类型问题）
if (nrow(outTab) > 0) {
  write.csv(outTab, "outTab.csv", row.names = FALSE)
  outTab <- read.csv("outTab.csv", header = TRUE)
  
  # 对结果进行分组标记
  outTab$group <- ifelse(outTab$pvalue > 0.05, "a",
                         ifelse(outTab$z > 0, "b", "c"))
  
  # 绘制Cox回归结果森林图
  p_cox <- ggplot(outTab, aes(x = gene, y = z, color = group)) + 
    geom_point(stat = 'identity', size = 3) +
    scale_color_manual(values = c("gray", "blue", "lightblue")) +
    geom_segment(aes(y = 0, x = gene, yend = z, xend = gene), 
                 color = "grey", linetype = "dashed") +
    scale_y_continuous(limits = c(min(c(-5, outTab$z)), max(c(5, outTab$z)))) +  # 自适应y轴范围
    theme_base() +
    theme(legend.position = "", axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip() +
    geom_hline(yintercept = c(2, -2), linetype = "dashed", colour = "grey", size = 1) +
    geom_hline(yintercept = c(3, -3), linetype = "solid") +
    xlab("") + ylab("Z score") +
    labs(title = "LUAD")
  
  print(p_cox)
} else {
  warning("No Cox regression results to plot")
}

save_plot("Fig/Fig6d_0.pdf", p_cox,base_height = 20,base_width = 10)

####################################################################################
# 生存曲线（KM曲线）分析部分

if (nrow(rt) > 0 && ncol(rt) >= 3) {
  # 复制数据用于生存曲线分析
  data_sur <- rt
  
  # 对每个细胞亚群的表达量进行高低分组（以中位数为界）
  for (j in 3:ncol(data_sur)) {
    # 仅对有足够数据的列进行分组
    if (sum(!is.na(data_sur[, j])) >= 2) {
      data_sur[, j] <- ifelse(as.numeric(data_sur[, j]) > median(data_sur[, j], na.rm = TRUE),
                              "high", "low")
    } else {
      data_sur[, j] <- NA
      warning(paste("Column", colnames(data_sur)[j], "has insufficient data for grouping"))
    }
  }
  
  # 创建空列表存储生存曲线
  plot <- list()
  valid_clusters <- 0  # 计数有效亚群
  
  # 循环绘制每个细胞亚群的KM曲线（替换原循环部分）
  for (j in 1:(ncol(data_sur) - 2)) {  # 排除前2列（survival和status）
    cluster_name <- colnames(data_sur)[j + 2]
    # 提取当前亚群的生存数据
    data_sur_plot <- data_sur[, c(1, 2, j + 2)]
    colnames(data_sur_plot) <- c("survival", "status", "group")
    
    # 过滤缺失值
    data_sur_plot <- data_sur_plot[!is.na(data_sur_plot$group) & 
                                     !is.na(data_sur_plot$survival) & 
                                     !is.na(data_sur_plot$status), ]
    
    # 确保两组都有样本才绘图
    if (nrow(data_sur_plot) >= 10 && length(unique(data_sur_plot$group)) == 2) {
      valid_clusters <- valid_clusters + 1
      # 绘制KM曲线，明确指定返回的绘图元素
      km <- ggsurvplot(
        survfit(Surv(survival, status) ~ group, data = data_sur_plot),
        pval = TRUE,
        legend = "right",
        legend.labs = c("Low", "High"),
        xlab = "Time",
        ylab = "Survival Probability",
        main = cluster_name,
        return = TRUE  # 关键：返回可处理的对象
      )
      # 提取曲线的绘图对象（grob），而非列表
      plot[[valid_clusters]] <- km$plot  # 直接提取plot组件
    }
  }
  
  # 组合并保存生存曲线
  if (valid_clusters > 0) {
    p.all <- cowplot::plot_grid(plotlist = plot, ncol = 2)
    # 创建Fig文件夹（若不存在）
    if (!dir.exists("Fig")) dir.create("Fig")
    save_plot("Fig/Fig6d.pdf", p.all, base_height = 3 * ceiling(valid_clusters / 2), base_width = 10)
    print("生存曲线已保存至Fig/Fig6d.pdf")
  } else {
    warning("No valid clusters for KM curve analysis")
  }
} else {
  warning("Insufficient data for survival curve analysis")
}


####################################################################################
# 编程猫中绘制图的代码 --------------------------------------------------------------
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
save_plot("Fig/Fig6d_1.pdf", p.all,base_height = 20,base_width = 10)

