
## http://tisch.compbio.cn/home/, 190 数据集 6,297,320 个细胞

# 1. 加载必要的包
pacman::p_load(
  Seurat,         # 单细胞数据处理
  dplyr,          # 数据处理
  fgsea,          # GSEA分析
  here,           # 文件路径处理
  data.table,    # 数据操作
  clusterProfiler # GSEA相关函数
)
# 2. 设置工作目录
here()

# 3. 读取 Hallmark 基因集（GMT 文件）
