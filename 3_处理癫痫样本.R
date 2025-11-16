# 加载必需包（仅保留核心依赖）
pacman::p_load(Seurat, dplyr, Matrix)

# 读取并初始化Seurat对象（核心功能：读取数据+初步过滤）
read_and_init_seurat <- function(sample_name, prefix, data_dir) {
  # 构建文件路径
  barcode_path <- file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz"))
  feature_path <- file.path(data_dir, paste0(prefix, "_features.tsv.gz"))
  matrix_path <- file.path(data_dir, paste0(prefix, "_matrix.mtx.gz"))
  
  # 检查文件存在性，不完整则跳过
  if (!all(file.exists(barcode_path, feature_path, matrix_path))) {
    return(NULL)
  }
  
  # 读取数据并处理
  counts <- readMM(gzfile(matrix_path))
  barcodes <- readLines(gzfile(barcode_path))
  features <- read.table(gzfile(feature_path), header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- features[, 2]
  colnames(counts) <- barcodes
  counts <- as(counts, "CsparseMatrix")
  counts <- counts[!duplicated(rownames(counts)), ]  # 去重基因
  
  # 创建Seurat对象（初步过滤：min.cells=3，min.features=200）
  sc_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # 计算线粒体基因比例（质控必需）
  sc_obj <- PercentageFeatureSet(sc_obj, pattern = "^MT-", col.name = "percent.mt")
  return(sc_obj)
}

# 质控函数（仅保留筛选高质量细胞，删除所有可视化）
qc_and_visualize <- function(sc_obj, sample_name, output_dir) {
  # 筛选高质量细胞（核心阈值不变）
  high_quality <- sc_obj@meta.data$nFeature_RNA > 200 & 
    sc_obj@meta.data$nFeature_RNA < 5000 & 
    sc_obj@meta.data$percent.mt < 10
  sc_hq <- subset(sc_obj, cells = rownames(sc_obj@meta.data[high_quality, ]))
  
  # 保存质控后对象（仅保留核心文件）
  save_path <- file.path(output_dir, sample_name)
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sc_hq, file.path(save_path, paste0(sample_name, "_hq.rds")))
  
  return(sc_hq)
}

# 主流程（批量处理+合并）
main <- function() {
  # 设置路径
  data_dir <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/1_raw_data/GSE302285_epilepsy"
  output_dir <- "/Users/lijiangbo/1_Projects/Epi1106/5_scRNA/2_processed_data/1113"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 样本信息
  samples <- list(
    "P010_frontal" = list(prefix = "GSM9101266_P010_frontal"),
    "P010_parietal" = list(prefix = "GSM9101267_P010_parietal"),
    "P018" = list(prefix = "GSM9101268_P018"),
    "P020" = list(prefix = "GSM9101269_P020"),
    "P024" = list(prefix = "GSM9101270_P024"),
    "P025" = list(prefix = "GSM9101271_P025")
  )
  
  # 批量处理样本（质控后保存）
  processed_objects <- list()
  for (sample_name in names(samples)) {
    prefix <- samples[[sample_name]]$prefix
    sc_obj <- read_and_init_seurat(sample_name, prefix, data_dir)
    if (!is.null(sc_obj)) {
      sc_hq <- qc_and_visualize(sc_obj, sample_name, output_dir)
      processed_objects[[sample_name]] <- sc_hq
    }
  }
  
  # 合并所有质控后样本
  if (length(processed_objects) > 0) {
    combined_obj <- merge(
      x = processed_objects[[1]],
      y = processed_objects[-1],
      add.cell.ids = names(processed_objects),
      project = "GSE302285_epilepsy6_combined"
    )
    # 保存合并后对象
    saveRDS(combined_obj, file.path(output_dir, "combined_epilepsy6_hq.rds"))
  }
}

# 执行
main()