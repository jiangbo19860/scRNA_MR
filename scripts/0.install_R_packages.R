## packages needed 
## author: 科研猫

## packages to install and manage R packages
install.packages("pacman")
install.packages("BiocManager")
install.packages("devtools")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  


## packges used in SC analysis
pacman::p_load(Seurat)
pacman::p_load(ggplot2)
pacman::p_load(tidyverse)
pacman::p_load(ggpubr)
pacman::p_load(patchwork)
pacman::p_load(GSVA)
pacman::p_load(clusterProfiler)
pacman::p_load(limma)
pacman::p_load(stringr)
pacman::p_load(pheatmap)
pacman::p_load(cowplot)
pacman::p_load(ggrepel)
pacman::p_load(survival)
pacman::p_load(ggthemes)
pacman::p_load(data.table)
pacman::p_load(survminer)

# check installation
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(pheatmap)
library(cowplot)
library(ggrepel)
library(survival)
library(ggthemes)
library(data.table)
library(survminer)
