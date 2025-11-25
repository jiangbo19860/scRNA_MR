# 设置文件路径
file1 <- "F:\\gongbing\\1.data\\ebi-a-GCST90025990_calcium_levels.vcf.gz"
file2 <- "F:\\gongbing\\1.data\\GCST90020053_buildGRCh37.tsv"

# 读取TSV文件
library(data.table)
data1 <- fread(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data2 <- fread(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(data1,130)
head(data2)

# ========== 步骤1：统一data1的结构（解决V1不存在的问题） ==========
# 如果data1是字符向量，转为单列data.table并命名为V1
if (is.vector(data1)) {
  data1 <- data.table(V1 = data1)
} else if (is.data.table(data1) && !"V1" %in% names(data1)) {  # else与上一行的}紧接
  setnames(data1, names(data1)[1], "V1")
}

# ========== 步骤2：过滤VCF header行（以##开头） ==========
data1_clean <- data1[!grepl("^##", V1)]
head(data1_clean)

# ========== 步骤3：拆分VCF数据行（按制表符分割列） ==========
data1_clean[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DATA") :=
              tstrsplit(V1, "\t", fixed = TRUE, fill = TRUE)]

# ========== 步骤4：拆分DATA列（ES:SE:LP:AF:ID） ==========
data1_clean[, c("ES", "SE", "LP", "AF", "ID_DATA") :=
              tstrsplit(DATA, ":", fixed = TRUE, fill = TRUE)]

# ========== 步骤5：数据类型转换与字段映射 ==========
data1_clean[, `:=`(
  variant_id = fifelse(ID == "." | is.na(ID), ID_DATA, ID),
  chromosome = as.integer(CHROM),
  base_pair_location = as.integer(POS),
  effect_allele = ALT,
  other_allele = REF,
  effect_allele_frequency = as.numeric(AF),
  beta = as.numeric(ES),
  standard_error = as.numeric(SE),
  LP_num = as.numeric(LP)
)]

# ========== 步骤6：计算p值（处理特殊值） ==========
data1_clean[, p_value := fifelse(
  is.infinite(LP_num) | LP_num > 300,
  0,
  10^(-LP_num)
)]
data1_clean[is.na(LP_num), p_value := NA_real_]

# ========== 步骤7：构建与data2一致的格式 ==========
data1_converted <- data1_clean[, .(
  variant_id,
  chromosome,
  base_pair_location,
  effect_allele,
  other_allele,
  effect_allele_frequency,
  beta,
  standard_error,
  p_value,
  Direction = NA_character_,
  HetISq = NA_real_,
  HetChiSq = NA_real_,
  HetPVal = NA_real_
)]

# 按染色体和位置排序
setorder(data1_converted, chromosome, base_pair_location)

# 验证结果
head(data1_converted)

# 保存修改后的数据为新的文件
output_file <- "/Volumes/AI_blue/gongbing/1.data/ebi-a-GCST90025990_calcium_levels.tsv"  # 更新文件保存路径
fwrite(data1, output_file, sep = "\t", quote = FALSE)

cat("数据已保存至：", output_file, "\n")
