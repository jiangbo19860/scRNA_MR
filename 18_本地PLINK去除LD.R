# 对 GWAS（全基因组关联研究）中的暴露数据（exposure data）进行连锁不平衡（LD）聚类（clumping），从数据中筛选出独立的、不相关的 SNP（单核苷酸多态性），为后续的孟德尔随机化（MR）分析或其他遗传关联分析做准备。
# 要先把PLINK软件下载到本地（https://www.cog-genomics.org/plink/1.9/），并配置好环境变量，确保R可以调用PLINK命令行工具。
# 还要下载LD 参考数据集（local_ld_panel）：指定本地的 1000 Genomes 项目参考数据（欧洲人群 EUR），包含.bed/.bim/.fam格式文件，用于计算 SNP 之间的连锁不平衡（LD）程度。

rm(list = ls())  # 清空工作空间
pacman::p_load(
  here,
  TwoSampleMR,
  ieugwasr,
  dplyr,
  data.table
)
here()

# 1. 设置工作目录（与数据存放路径一致）
setwd(here("4_references/Stroke_MR+机器学习/37去除连锁不平衡性"))

# 3. 读取暴露数据
exposure_data <- read_exposure_data(
  filename = "exposure_data_filtered.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  pval_col = "pval.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  samplesize_col = "samplesize.exposure",
  chr_col = "chr.exposure",
  pos_col = "pos.exposure",
  clump = FALSE  # 关闭自动聚类，手动处理
)

# 4. 配置本地资源（关键步骤）
plink_abs_path <- "/Users/lijiangbo/bin/plink"
file.exists(plink_abs_path)  # 返回 TRUE 表示文件存在
file.access(plink_abs_path, mode = 1) == 0  # 返回 TRUE 表示有执行权限，mode = 0：检查文件是否存在； mode = 1：检查是否有执行权限， mode = 2：检查是否有写入权限； mode = 4：检查是否有读权限。
# 调用 PLINK 并获取版本信息
plink_version <- tryCatch(
  system2(plink_abs_path, args = "--version", stdout = TRUE, stderr = TRUE),
  error = function(e) e  # 捕获错误信息
)
print(plink_version)

# 4.2 本地LD参考数据路径（PLINK格式，需提供.bed/.bim/.fam文件的前缀）
# 确保该路径下存在：1000G_EUR.bed、1000G_EUR.bim、1000G_EUR.fam三个文件
local_ld_panel <- "/Users/lijiangbo/scRNA_MR/1_data/1kg.v3/EUR"  # 无需加文件后缀

# 5. 执行LD聚类（使用本地数据和PLINK，核心修正）
exposure_dat_clumped <- clump_data(
  dat = exposure_data,
  clump_kb = 10000,   # 聚类窗口：10000kb（10Mb）即只对 10Mb 范围内的 SNP 进行 LD 过滤。
  clump_r2 = 0.001,        # LD阈值：r² > 0.001的SNP会被剔除
  bfile = local_ld_panel,  # 本地参考数据集（核心参数，替代API）
  plink_bin = plink_abs_path  # 指定本地 PLINK 软件，用于执行 LD 计算。
)

# 6. 查看聚类结果（验证是否成功）
cat("原始SNP数量：", nrow(exposure_data), "\n")
cat("聚类后保留的SNP数量：", nrow(exposure_dat_clumped), "\n")

# 7. 保存处理后的结果
write.csv(
  exposure_dat_clumped,
  file = here("3_outputs/exposure_clumped.csv"),  # 使用 here() 定位到 3_outputs 文件夹
  row.names = FALSE
)
