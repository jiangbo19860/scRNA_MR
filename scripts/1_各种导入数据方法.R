rm(list = ls())
pacman::p_load(Seurat, dplyr, data.table, ggplot2, cowplot, patchwork, ggrepel, hdf5r)
# 设置工作目录
here::i_am("scRNAseq.Rproj")

###### 第一种：读取10X Genomics的标准格式数据 #####
pbmc.data <- Read10X(data.dir = "scRNA_read/filtered_feature_bc_matrix/")
# 创建Seurat对象，这是单细胞数据分析的核心数据结构
# 作用：将原始表达矩阵整合为结构化对象，同时进行初步的数据过滤
pbmc <- CreateSeuratObject(
  counts = pbmc.data, # 输入的基因表达计数矩阵（行=基因，列=细胞，值=UMI计数）
  # 通常来自10X Genomics的filtered_feature_bc_matrix文件夹
  project = "pbmc500", # 项目名称（自定义，用于标识数据集，例如"pbmc500"代表500个PBMC细胞）
  # 会存储在对象的元数据中，方便多数据集对比
  min.cells = 3, # 基因过滤阈值：仅保留在至少3个细胞中表达的基因
  # 目的：剔除仅在1-2个细胞中表达的基因（可能是技术噪声或测序错误）
  min.features = 200 # 细胞过滤阈值：仅保留检测到至少200个基因的细胞
  # 目的：剔除低质量细胞（如死细胞，因细胞膜破裂导致基因检测数少）
)

# 查看创建的Seurat对象信息（运行后会显示对象的核心属性）
pbmc

#### 第二种：读取10X Genomics的h5格式数据，也可以得到一样的数据 #####
# pbmc.data1是变量名+数据类型标识，是原始的稀疏表达矩阵（仅包含基因表达数据）
pbmc.data1 <- Read10X_h5(filename = "scRNA_read/500_PBMC_3p_LT_Chromium_X_filtered_feature_bc_matrix.h5")
pbmc1 <- CreateSeuratObject(
  counts = pbmc.data1,
  project = "pbmc_h5",
  min.cells = 3,
  min.features = 200
)

###### 第三种读取方式：依次读取每个子文件，使用 readMM() 读取压缩文件 #####
library(Matrix)

# 读取表达矩阵（注意使用 gzfile() 处理压缩文件）
counts <- readMM(gzfile("scRNA_read/matrix.mtx.gz"))

# 读取细胞条形码
barcodes <- readLines(gzfile("scRNA_read/barcodes.tsv.gz"))

# 读取基因信息
features <- read.table(gzfile("scRNA_read/features.tsv.gz"),
  header = FALSE,
  stringsAsFactors = FALSE
)

# 为基因表达矩阵添加行名和列名，使其成为一个“带注释的”表达矩阵。
rownames(counts) <- features[, 2] # 1表示基因ID（ENSG开头的）, 2表示基因名。
colnames(counts) <- barcodes # 细胞条形码

# 查看矩阵结构
dim(counts) # 输出：基因数 x 细胞数
counts[1:5, 1:5]
counts <- as(counts, "CsparseMatrix") # 推荐写法
class(counts)


###### 4. 从GEO数据库中读取 #####
# 压缩的csv文件
data <- read.table(gzfile("scRNA_read/GSE158631_count.csv.gz"), sep = ",", header = TRUE)
data <- data |>
  column_to_rownames("X") |> # 将第一列设置为行名
  as.matrix() |> # 转换为矩阵
  as("sparseMatrix") # 转换为稀疏矩阵
