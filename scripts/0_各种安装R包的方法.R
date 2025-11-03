
# 1. 确保 BiocManager 是最新版本
install.packages("BiocManager")

# 2. 升级 Bioconductor 到 3.21（与 R 4.5 兼容）
BiocManager::install(version = "3.21")


# 1. 安装下载到本地的R包-impute -----------------------------------------------------------
# 安装本地的impute包，repos = NULL表示不使用在线仓库，type = "source"表示安装源包
install.packages("/Users/lijiangbo/scRNAseq/R包/impute_1.82.0.tgz", repos = NULL, type = "source")

install.packages("/Users/lijiangbo/scRNAseq/R包/BiocVersion_3.21.1.tgz", repos = NULL, type = "source")  # https://bioconductor.org/packages/release/bioc/html/preprocessCore.html

if (!require("remotes")) {
  install.packages("remotes")  # 从 CRAN 安装 remotes
}
remotes::install_local("/Users/lijiangbo/scRNAseq/R包/presto-master.zip")

remotes::install_local("/Users/lijiangbo/scRNAseq/R包/RSQLite-main.zip")

install.packages("/Users/lijiangbo/scRNAseq/R包/RSQLite_2.4.1.tgz", repos = NULL, type = "binary")  # 安装RSQLite包

remotes::install_github("satijalab/azimuth")

# 2. 安装Bioconductor包--https://bioconductor.org/packages/release/bioc/html/preprocessCore.html -----------------------------------------------------------
# 从Bioconductor安装
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::available()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocVersion")
library(BiocVersion)

# 3. 安装Github中的R包-WGCNA -----------------------------------------------------------
# 1. 确保安装了 devtools 或 remotes 包（用于从 GitHub 安装）
if (!require("remotes")) {
  install.packages("remotes")  # 从 CRAN 安装 remotes
}

# 2. 从 GitHub 安装 presto
remotes::install_github("immunogenomics/presto")
install.packages('RSQLite')

# install BiocManager
install.packages("BiocManager")

# install Bioconductor core packages
BiocManager::install()

# install devtools
BiocManager::install("devtools")

# install additional packages
BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))

library(WGCNA)

# 安装包（如果未安装）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "GO.db"))

# 加载包
library(AnnotationDbi)
library(GO.db)

# 安装和加载必要的包
BiocManager::install("Rgraphviz")
library(Rgraphviz)
# 安装CRAN上的RSQLite包
install.packages("RSQLite")


BiocManager::install("BiocCheck")



# 确保BiocManager已安装且更新
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 更换为稳定的镜像源（如中科大）
options(
  repos = c(
    CRAN = "https://mirrors.ustc.edu.cn/CRAN/",
    BiocPackages = "https://mirrors.ustc.edu.cn/bioc/"
  )
)

# 安装preprocessCore包
BiocManager::install("preprocessCore")

# 再次尝试加载WGCNA
library(WGCNA)



# 4. renv::snapshot()时要求安装的R包 ---------------------------------------------
# R包安装脚本 - 用于批量安装所需依赖
# 注意：运行前确保已安装BiocManager包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 检查renv包是否安装
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# 定义需要安装的CRAN包列表
cran_packages <- c(
  "cowplot", "devtools", "dlm", "forcats", "future", "ggpubr", 
  "ggrepel", "ggridges", "ggsci", "ggthemes", "here", "hdf5r", 
  "MCMCpack", "mixtools", "pacman", "parallelDist", "patchwork", 
  "pheatmap", "readr", "remotes", "reshape2", "showtext", "survminer", 
  "sysfonts", "tidyverse", "umap", "usethis"
)

# 定义需要安装的Bioconductor包列表
bioc_packages <- c(
  "AnnoProbe", "BiocParallel", "celldex", "clusterProfiler", 
  "enrichplot", "fgsea", "GSVA", "limma", "loomR", "monocle", 
  "monocle3", "org.Hs.eg.db", "presto", "scater", "SCopeLoomR", 
  "scran", "sctransform", "SingleCellExperiment", "SingleR"
)

# 定义需要从GitHub安装的包列表
github_packages <- c(
  "satijalab/seurat",         # Seurat
  "satijalab/azimuth",        # Azimuth (GitHub源)
  "chris-mcginnis-ucsf/infercnv",  # infercnv
  "Vivianstats/copykat",      # copykat
  "saezlab/dorothea",        # dorothea
  "saezlab/nichenetr",       # nichenetr
  "chris-mcginnis-ucsf/IOBR", # IOBR
  "chris-mcginnis-ucsf/harmony", # harmony
  "sqjin/CellChat"            # CellChat
)

# 安装CRAN包
install_cran_packages <- function(packages) {
  cat("正在安装CRAN包...\n")
  for (pkg in packages) {
    cat(paste0("安装 ", pkg, "...\n"))
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        cat(paste0("✓ ", pkg, " 安装成功\n"))
      }, error = function(e) {
        cat(paste0("✗ ", pkg, " 安装失败: ", e$message, "\n"))
      })
    } else {
      cat(paste0("✓ ", pkg, " 已安装\n"))
    }
  }
}

# 安装Bioconductor包
install_bioc_packages <- function(packages) {
  cat("正在安装Bioconductor包...\n")
  for (pkg in packages) {
    cat(paste0("安装 ", pkg, "...\n"))
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        cat(paste0("✓ ", pkg, " 安装成功\n"))
      }, error = function(e) {
        cat(paste0("✗ ", pkg, " 安装失败: ", e$message, "\n"))
      })
    } else {
      cat(paste0("✓ ", pkg, " 已安装\n"))
    }
  }
}

# 安装GitHub包
install_github_packages <- function(packages) {
  cat("正在安装GitHub包...\n")
  for (pkg in packages) {
    pkg_name <- tail(strsplit(pkg, "/")[[1]], 1)
    cat(paste0("安装 ", pkg_name, " 从 ", pkg, "...\n"))
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      tryCatch({
        devtools::install_github(pkg, dependencies = TRUE)
        cat(paste0("✓ ", pkg_name, " 安装成功\n"))
      }, error = function(e) {
        cat(paste0("✗ ", pkg_name, " 安装失败: ", e$message, "\n"))
      })
    } else {
      cat(paste0("✓ ", pkg_name, " 已安装\n"))
    }
  }
}

# 执行安装
install_cran_packages(cran_packages)
install_bioc_packages(bioc_packages)
install_github_packages(github_packages)

# 验证所有包是否安装成功
check_packages <- function(packages) {
  missing <- c()
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  return(missing)
}

# 检查CRAN包
missing_cran <- check_packages(cran_packages)
# 检查Bioconductor包
missing_bioc <- check_packages(bioc_packages)
# 检查GitHub包(只检查包名部分)
github_pkg_names <- sapply(github_packages, function(x) tail(strsplit(x, "/")[[1]], 1))
missing_github <- check_packages(github_pkg_names)

# 输出安装结果
cat("\n===== 安装结果汇总 =====\n")
if (length(missing_cran) == 0 && length(missing_bioc) == 0 && length(missing_github) == 0) {
  cat("✓ 所有包安装成功!\n")
} else {
  cat("✗ 以下包安装失败:\n")
  if (length(missing_cran) > 0) cat("  CRAN包: ", paste(missing_cran, collapse = ", "), "\n")
  if (length(missing_bioc) > 0) cat("  Bioconductor包: ", paste(missing_bioc, collapse = ", "), "\n")
  if (length(missing_github) > 0) cat("  GitHub包: ", paste(missing_github, collapse = ", "), "\n")
}

# 更新renv快照(可选)
if (interactive()) {
  cat("\n是否要更新renv快照? (y/n): ")
  answer <- readline()
  if (tolower(answer) == "y") {
    tryCatch({
      renv::snapshot()
      cat("✓ renv快照已更新\n")
    }, error = function(e) {
      cat("✗ 更新renv快照失败: ", e$message, "\n")
    })
  }
}




# 5. 检查包是否已经安装 ---------------------------------------------------------------
# 检查单个包（返回逻辑值TRUE/FALSE）
require("AnnoProbe", quietly = TRUE)  # 不显示加载信息

# 批量检查多个包
required_packages <- c(
  "AnnoProbe", "AnnotationDbi", "Azimuth", "Biobase", "BiocParallel", 
  "BiocVersion", "CellChat", "celldex", "clusterProfiler", "copykat", 
  "doParallel", "dorothea", "enrichplot", "fgsea", "foreach", "future", 
  "GO.db", "GSVA", "harmony", "here", "infercnv", "IOBR", "limma", 
  "loomR", "monocle", "monocle3", "nichenetr", "org.Hs.eg.db", "pheatmap", 
  "presto", "readr", "Rgraphviz", "scater", "SCopeLoomR", "scran", 
  "sctransform", "Seurat", "SingleCellExperiment", "SingleR", "tidyverse", 
  "umap", "viridis", "WGCNA"
)

# 检查并输出每个包的安装状态
package_status <- sapply(required_packages, function(pkg) {
  require(pkg, character.only = TRUE, quietly = TRUE)
})

# 显示未安装的包
cat("未安装的包：\n")
print(required_packages[!package_status])



# 0720安装 ------------------------------------------------------------------
# 定义需要安装的包
missing_packages <- c(
  "AnnoProbe", "Azimuth", "BiocParallel", "BiocVersion", 
  "CellChat", "clusterProfiler", "copykat", "dorothea", 
  "enrichplot", "fgsea", "GO.db", "GSVA", 
  "infercnv", "IOBR", "limma", "loomR", 
  "monocle", "monocle3", "nichenetr", "org.Hs.eg.db", 
  "presto", "Rgraphviz", "SCopeLoomR", "Seurat", 
  "viridis", "WGCNA"
)

# 安装BiocManager（如果未安装）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 安装remotes（用于从GitHub安装包）
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 检查并更新Bioconductor版本
BiocManager::install(version = "3.21")  # 根据你的R版本选择合适的Bioconductor版本

# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 创建安装进度记录
install_log <- data.frame(
  package = missing_packages,
  source = NA,
  status = NA,
  message = NA,
  stringsAsFactors = FALSE
)

# 按来源分类包（根据实际情况调整）
bioc_packages <- c(
  "AnnoProbe", "BiocParallel", "BiocVersion", "CellChat", 
  "clusterProfiler", "dorothea", "enrichplot", "fgsea", 
  "GO.db", "GSVA", "infercnv", "IOBR", "limma", 
  "loomR", "monocle", "monocle3", "nichenetr", "org.Hs.eg.db", 
  "Rgraphviz", "SCopeLoomR", "Seurat"
)

cran_packages <- c("viridis", "presto")

github_packages <- c(
  "Azimuth" = "satijalab/azimuth",  # GitHub格式：用户名/仓库名
  "WGCNA" = "cran/WGCNA"  # 某些包需要特殊处理
)

# 安装Bioconductor包
for (pkg in bioc_packages) {
  if (pkg %in% missing_packages) {
    cat(paste0("正在安装Bioconductor包: ", pkg, "\n"))
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      install_log[install_log$package == pkg, c("source", "status")] <- c("Bioconductor", "成功")
    }, error = function(e) {
      msg <- paste("错误:", e$message)
      install_log[install_log$package == pkg, c("source", "status", "message")] <- c("Bioconductor", "失败", msg)
      cat(paste0("  ❌ 安装失败: ", pkg, " (", msg, ")\n"))
    })
  }
}

# 安装CRAN包
for (pkg in cran_packages) {
  if (pkg %in% missing_packages) {
    cat(paste0("正在安装CRAN包: ", pkg, "\n"))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
      install_log[install_log$package == pkg, c("source", "status")] <- c("CRAN", "成功")
    }, error = function(e) {
      msg <- paste("错误:", e$message)
      install_log[install_log$package == pkg, c("source", "status", "message")] <- c("CRAN", "失败", msg)
      cat(paste0("  ❌ 安装失败: ", pkg, " (", msg, ")\n"))
    })
  }
}

# 安装GitHub包
for (pkg in names(github_packages)) {
  if (pkg %in% missing_packages) {
    repo <- github_packages[pkg]
    cat(paste0("正在安装GitHub包: ", pkg, " (来源: ", repo, ")\n"))
    tryCatch({
      remotes::install_github(repo, dependencies = TRUE)
      install_log[install_log$package == pkg, c("source", "status")] <- c("GitHub", "成功")
    }, error = function(e) {
      msg <- paste("错误:", e$message)
      install_log[install_log$package == pkg, c("source", "status", "message")] <- c("GitHub", "失败", msg)
      cat(paste0("  ❌ 安装失败: ", pkg, " (", msg, ")\n"))
    })
  }
}

# 显示安装结果
cat("\n===== 安装结果汇总 =====\n")
print(install_log)

# 检查未成功安装的包
failed_packages <- install_log$package[install_log$status == "失败"]
if (length(failed_packages) > 0) {
  cat("\n以下包安装失败，请手动处理:\n")
  print(failed_packages)
} else {
  cat("\n所有包均安装成功！\n")
}

