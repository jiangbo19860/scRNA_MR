library(data.table)
library(dplyr)

# --------------------- 1. 仅使用 CPASSOC 数据（去掉MTAG相关逻辑） ---------------------
# 读取CPASSOC数据（SNP列 + p_SHet列）
dat_CPASSOC <- data.table::fread("./CPASSOC/cpassoc_sig.txt")[, c(1, 8)]  # 仅保留SNP和p_SHet
cat("CPASSOC 原始数据维度：", dim(dat_CPASSOC), "\n")
head(dat_CPASSOC)

# 筛选 p_SHet < 5e-8 的显著SNP
dat_CPASSOC <- subset(dat_CPASSOC, p_SHet < 5e-8)
cat("CPASSOC 筛选后数据维度：", dim(dat_CPASSOC), "\n")
head(dat_CPASSOC)

# 直接使用筛选后的CPASSOC数据作为分析数据（无MTAG合并）
data <- dat_CPASSOC
# 关键：添加后续筛选需要的列（如果CPASSOC没有单独trait的p值，用p_SHet填充）
data[, c("trait1_pval", "trait2_pval") := p_SHet]
# 占位染色体/位置列（后续HESS匹配需补充实际信息，无则保留NA）
if (!"trait1_CHR" %in% colnames(data)) data[, trait1_CHR := NA_integer_]
if (!"trait1_BP" %in% colnames(data)) data[, trait1_BP := NA_integer_]

# --------------------- 2. LD聚类分析（修复参数传递：去掉 data = ） ---------------------
sig_ind_snp <- clump_data_local_Online(
  data,  # 直接传入数据，无需指定 data = 参数名
  snp_col = "SNP",
  pval_col = "p_SHet",
  clump_pval = 5e-8,
  clump_kb = 500,
  clump_r2 = 0.3,
  bfile_1000G = "./1.data/1000G/EUR"
)
cat("LD聚类后剩余SNP数量：", nrow(sig_ind_snp), "\n")
head(sig_ind_snp)

# --------------------- 3. 筛选与novel_SNPs定义 ---------------------
# 无需trait1/trait2 pval筛选（直接用LD聚类结果）
sig_pre_SNP <- sig_ind_snp

# 进一步筛选高LD区域的SNP（clump_LD_R2同样修复参数传递）
r2_SNP <- clump_LD_R2(
  sig_pre_SNP,  # 直接传入数据，无需参数名
  snp_col = "SNP",
  clump_kb = 1000,
  clump_r2 = 0.3,
  bfile_1000G = "./1.data/1000G/EUR"
)

# 排除高LD的SNP
dat3 <- subset(sig_ind_snp, !SNP %in% r2_SNP)
cat("排除高LD后SNP数量：", nrow(dat3), "\n")
# novel_SNPs 维度： 58 6
# 排除排除高LD后SNP数量： 58

# 定义novel_SNPs（直接用dat3，因无MTAG相关筛选）
novel_SNPs <- dat3
cat("novel_SNPs 维度：", dim(novel_SNPs), "\n")
head(novel_SNPs)

# --------------------- 4. （可选）与HESS结果合并 ---------------------
if (file.exists("./HESS/HESS/step3.txt")) {
  hess_dat <- data.table::fread("./HESS/HESS/step3.txt")
  hess_dat$fdr <- p.adjust(hess_dat$p, method = "bonferroni")
  hess_dat <- subset(hess_dat, fdr < 0.05)
  cat("HESS 显著区域数量：", nrow(hess_dat), "\n")

  novel_SNPs_merge <- data.table()
  if (nrow(hess_dat) > 0 && nrow(novel_SNPs) > 0) {
    # 检查是否有有效染色体/位置信息
    if (all(is.na(novel_SNPs$trait1_CHR)) || all(is.na(novel_SNPs$trait1_BP))) {
      warning("novel_SNPs 缺少染色体/位置信息，无法与HESS区域匹配！")
    } else {
      for (j in 1:nrow(hess_dat)) {
        hess_data <- hess_dat[j, ]
        novel_SNPs_new <- novel_SNPs %>%
          dplyr::filter(
            trait1_CHR == hess_data$chr,
            trait1_BP > hess_data$start,
            trait1_BP < hess_data$end
          )
        novel_SNPs_merge <- rbind(novel_SNPs_merge, novel_SNPs_new)
      }
      cat("HESS区域匹配后的SNP数量：", nrow(novel_SNPs_merge), "\n")
      head(novel_SNPs_merge)
    }
  }
} else {
  warning("未找到HESS文件，跳过HESS合并步骤")
  novel_SNPs_merge <- data.table()
}
# HESS 显著区域数量： 0
# 保存最终结果
fwrite(novel_SNPs, "./novel_SNPs_final.txt", sep = "\t")
if (nrow(novel_SNPs_merge) > 0) {
  fwrite(novel_SNPs_merge, "./novel_SNPs_HESS_merge.txt", sep = "\t")
}
