rm(list = ls())
# 1. 安装并加载所需包
pacman::p_load(
  here,        # 路径管理
  readr,       # 高效读取文件
  dplyr,       # 数据处理
  tidyr,       # 缺失值处理（提供drop_na函数）
  TwoSampleMR, # MR分析核心包
  ieugwasr     # 辅助GWAS数据处理
)
here()

# 1. 定义压缩文件路径
compressed_file <- here("1_data/finngen/finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz")

# 检查压缩文件是否存在
if (!file.exists(compressed_file)) {
  stop("压缩文件不存在，请检查路径：", compressed_file)
}

# 2. 定义解压后的临时文件名（无后缀）和目标文件名（带.tsv）
temp_uncompressed <- tools::file_path_sans_ext(compressed_file, compression = TRUE)  # 无后缀的解压文件名
uncompressed_file <- paste0(temp_uncompressed, ".tsv")  # 最终目标文件名（带.tsv）

# 3. 解压文件（得到无后缀的临时文件）
if (!file.exists(uncompressed_file)) {  # 仅当目标文件不存在时操作
  # 解压.gz文件（得到无后缀的文件）
  cat("正在解压文件...\n")
  R.utils::gunzip(
    filename = compressed_file,
    destname = temp_uncompressed,  # 先解压为无后缀文件
    overwrite = TRUE,
    remove = FALSE
  )

  # 4. 检查临时文件是否存在，然后添加.tsv后缀
  if (file.exists(temp_uncompressed)) {
    # 重命名文件，添加.tsv后缀
    file.rename(from = temp_uncompressed, to = uncompressed_file)
    cat("已为文件添加.tsv后缀：", uncompressed_file, "\n")
  } else {
    stop("解压失败，未生成临时文件：", temp_uncompressed)
  }
} else {
  cat("目标文件已存在（带.tsv后缀）：", uncompressed_file, "\n")
}

# 5. 验证最终文件
if (file.exists(uncompressed_file)) {
  cat("文件准备完成，开始读取列名...\n")
  col_names <- readr::read_tsv(uncompressed_file, n_max = 1, show_col_types = FALSE) %>% colnames()
  cat("文件列名：", paste(col_names, collapse = ", "), "\n")
} else {
  stop("最终文件不存在，请检查操作：", uncompressed_file)
}

# 4. 读取并整理数据
## 4.1 读取解压后的文件（FinnGen通常为TSV格式）
# 查看数据列名（根据实际列名调整，以下为常见列名示例）
col_names <- readr::read_tsv(uncompressed_file, n_max = 1) %>% colnames()
cat("数据列名：", paste(col_names, collapse = ", "), "\n")

# 检查文件是否存在
if (!file.exists(uncompressed_file)) {
  stop("解压后的文件不存在，请检查路径：", uncompressed_file)
}

# 读取数据并整理为MR所需格式
finngen_data <- readr::read_tsv(
  uncompressed_file,
  col_types = readr::cols(),  # 自动识别列类型
  show_col_types = FALSE     # 隐藏列类型提示
) %>%
  # 重命名为TwoSampleMR要求的标准列名（根据实际列名调整）
  dplyr::rename(
    SNP = rsids,               # SNP标识符（对应数据中的rsids列）
    beta = beta,               # 效应值（对应数据中的beta列）
    se = sebeta,               # 标准误（对应数据中的sebeta列）
    pval = pval,               # P值（对应数据中的pval列）
    effect_allele = alt,       # 效应等位基因（对应数据中的alt列）
    other_allele = ref,        # 非效应等位基因（对应数据中的ref列）
    eaf = af_alt,              # 效应等位基因频率（对应数据中的af_alt列）
    chr = `#chrom`,            # 染色体（对应数据中的#chrom列，注意#需要用`包裹）
    pos = pos                  # 位置（对应数据中的pos列）
  ) %>%
  # 保留MR分析必需的列
  dplyr::select(SNP, beta, se, pval, effect_allele, other_allele, eaf, chr, pos) %>%
  # 去除重复SNP（保留第一个出现的）
  dplyr::distinct(SNP, .keep_all = TRUE) %>%
  # 去除含有缺失值的行（使用tidyr包的drop_na，解决之前的错误）
  tidyr::drop_na()

# 查看处理后的数据量
cat("处理后的数据量：", nrow(finngen_data), "个SNP\n")
cat("前5行数据预览：\n")
print(head(finngen_data))

# 格式化为MR分析所需的暴露/结局数据（以结局数据为例，若为暴露则改为type = "exposure"）
outcome_dat <- TwoSampleMR::format_data(
  finngen_data,
  type = "outcome",           # 定义为结局数据（根据研究设计调整）
  phenotype_col = "glioma",   # 表型名称（胶质瘤）
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  chr_col = "chr",
  pos_col = "pos"
)

# 保存处理后的数据集
saveRDS(outcome_dat, here("3_outputs", "finngen_outcome_dat.rds"))
cat("数据处理完成，已保存至：", here("3_outputs", "finngen_outcome_dat.rds"), "\n")

head(outcome_dat)

# 6. LD聚类（如果是暴露数据，则进行下面的去除连锁不平衡的SNP），结局变量（outcome）的数据通常不需要进行 LD 聚类。
## 6.1 配置本地PLINK和参考数据集（需提前准备）
# plink_path <- "/Users/lijiangbo/bin/plink"  # PLINK可执行文件路径
# local_ld_panel <- here("1_data/1kg.v3/EUR")  # 1000G欧洲人群参考数据（前缀）
#
# ## 6.2 执行聚类
# exposure_clumped <- TwoSampleMR::clump_data(
#   exposure_dat,
#   clump_kb = 10000,
#   clump_r2 = 0.001,
#   bfile = local_ld_panel,
#   plink_bin = plink_path
# )
#
# cat("聚类后保留的SNP数量：", nrow(exposure_clumped), "\n")
#
# # 7. 保存处理后的数据集（用于后续MR分析）
# saveRDS(exposure_clumped, here("3_outputs/finngen_exposure_clumped.rds"))
# saveRDS(outcome_dat, here("3_outputs/finngen_outcome_dat.rds"))
# cat("处理完成，结果已保存到output文件夹\n")
