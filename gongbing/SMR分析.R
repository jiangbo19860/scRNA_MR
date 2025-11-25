# 1. 清空环境与加载必需包（整合示例代码的包需求）
rm(list = ls())
gc()
library(ggplot2)
library(dplyr)
library(stringr)
library(openxlsx)
library(purrr)
library(data.table)  # 新增：借鉴示例代码的数据处理包

# 2. 核心路径配置（融合示例代码的文件夹逻辑）
# 2.1 基础路径
adhd_gwas <- "/Volumes/AI_blue/gongbing/2.pre_data/ADHD_SMR.ma"       # 结局1：ADHD
epilepsy_gwas <- "/Volumes/AI_blue/gongbing/2.pre_data/epilepsy_SMR.ma" # 结局2：癫痫
gtex_eqtl_dir <- "/Volumes/AI_blue/gongbing/SMR/GTEx_brain/"           # GTEx 13脑区eQTL
smr_bin <- "/Users/lijiangbo/3_Resources/Bioinformatics/Tools/smr/smr" # SMR软件
ref_bfile <- "/Volumes/AI_blue/gongbing/MAGMA/g1000_eur/g1000_eur" # 参考基因组

# 2.2 结果输出路径（借鉴示例代码：创建子文件夹分类存放）
workpath <- "/Volumes/AI_blue/gongbing/3.SMR_results_GTEx13"
filename <- "SMR_result_files"  # 结果子文件夹名（示例代码逻辑）
result_dir <- file.path(workpath, filename)
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE) # 自动创建子文件夹
setwd(workpath)

# 2.3 定义结局列表（批量处理ADHD+癫痫，替代重复代码）
outcome_list <- list(
  list(name = "ADHD", path = adhd_gwas),
  list(name = "Epilepsy", path = epilepsy_gwas)
)

# 3. 自动识别13个脑区的eQTL文件（保留原逻辑，适配.lite格式）
esi_files <- list.files(
  path = gtex_eqtl_dir,
  pattern = "\\.lite\\.esi$",
  full.names = TRUE
)
if (length(esi_files) == 0) stop("❌ 未找到GTEx eQTL文件（需.lite.esi/.epi/.besd）")

# 提取脑区名称和eQTL文件前缀（如Brain_Amygdala.lite）
brain_regions <- map_dfr(esi_files, function(esi) {
  fname <- basename(esi)
  eqtl_prefix <- str_remove(fname, "\\.esi$")  # 得到Brain_Amygdala.lite
  region_name <- str_remove(eqtl_prefix, "\\.lite$")  # 得到Brain_Amygdala
  data.frame(
    region_name = region_name,
    eqtl_prefix = file.path(gtex_eqtl_dir, eqtl_prefix),  # eQTL文件前缀（带路径）
    stringsAsFactors = FALSE
  )
}) %>% distinct()

cat("✅ 自动识别到", nrow(brain_regions), "个GTEx脑区：\n")
print(brain_regions$region_name)

# 4. SMR运行函数
run_smr <- function(outcome_name, outcome_path, region_name, eqtl_prefix) {
  # 结果前缀（存子文件夹）
  res_prefix <- file.path(result_dir, paste0(outcome_name, "-", region_name))
  smr_result_file <- paste0(res_prefix, ".smr")
  
  # 1. 执行SMR命令（保留之前的优化参数）
  smr_cmd <- paste0(
    smr_bin,
    " --bfile '", ref_bfile, "'",          # 注意：仍需移至非iCloud路径！
    " --gwas-summary '", outcome_path, "'",
    " --beqtl-summary '", eqtl_prefix, "'",
    " --diff-freq 0.2",
    " --diff-freq-prop 0.90",
    " --maf 0.01",
    " --out '", res_prefix, "'",
    " --thread-num 6"
  )
  
  cat(paste0("\n=== 运行结局：", outcome_name, " | 脑区：", region_name, " ===\n"))
  cat("SMR命令：", smr_cmd, "\n")
  system(smr_cmd)  # 执行SMR
  
  # 2. 读取SMR结果并检查列名
  if (!file.exists(smr_result_file)) {
    cat("⚠️  未生成SMR结果文件\n")
    return(NULL)
  }
  smr_data <- read.delim(smr_result_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(smr_data) == 0) {
    cat("⚠️  SMR结果文件为空\n")
    return(NULL)
  }
 ##  colnames(smr_data)
  # [1] "probeID"    "ProbeChr"   "Gene"       "Probe_bp"   "topSNP"     "topSNP_chr" "topSNP_bp"  "A1"         "A2"        
  # [10] "Freq"       "b_GWAS"     "se_GWAS"    "p_GWAS"     "b_eQTL"     "se_eQTL"    "p_eQTL"     "b_SMR"      "se_SMR"    
  # [19] "p_SMR"      "p_HEIDI"    "nsnp_HEIDI"
  # 3. 动态识别P值列名
  # 3. 关键：提前定义所有列名变量（含p_HEIDI），并验证存在性
  # 基于断点确认的列名：p_SMR、b_SMR、se_SMR、p_HEIDI、topSNP
  p_col <- "p_SMR"               # P值列
  beta_col <- "b_SMR"           # 效应量列
  se_col <- "se_SMR"            # 标准误列
  p_heidi_col <- "p_HEIDI"      # HEIDI检验列（新增：提前定义）
  snp_col <- "topSNP"           # SNP列（SMR结果中SNP列名为topSNP）
  
  # 验证所有关键列是否存在（新增p_heidi_col检查）
  all_required_cols <- c(p_col, beta_col, se_col, p_heidi_col, snp_col, "Gene")
  missing_cols <- setdiff(all_required_cols, colnames(smr_data))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "❌ 以下关键列缺失：", paste(missing_cols, collapse = ", "), "\n",
      "SMR结果实际列名：", paste(colnames(smr_data), collapse = ", "), "\n",
      "请核对SMR结果列名与函数中定义的列名是否一致"
    ))
  }
  cat(paste0("✅ 确认所有列名：P值=", p_col, ", 效应量=", beta_col, ", 标准误=", se_col, ", HEIDI=", p_heidi_col, "\n"))
  
  # 4. 数据处理（所有列均用动态引用，避免硬编码错误）
  processed_data <- smr_data %>%
    mutate(
      # FDR校正（动态引用p_col）
      FDR = p.adjust(!!sym(p_col), method = "BH"),
      # OR及95%CI（动态引用beta_col/se_col）
      OR = exp(!!sym(beta_col)),
      Lower_CI = exp(!!sym(beta_col) - 1.96 * !!sym(se_col)),
      Upper_CI = exp(!!sym(beta_col) + 1.96 * !!sym(se_col)),
      OR_95CI = paste0(round(OR, 3), " (", round(Lower_CI, 3), "-", round(Upper_CI, 3), ")"),
      # 补充元信息
      Outcome = outcome_name,
      Brain_Region = region_name,
      EQTL_Prefix = basename(eqtl_prefix)
    ) %>%
    # 关键修改：p_HEIDI用动态引用!!sym(p_heidi_col)，而非直接写p_HEIDI
    select(
      Outcome, Brain_Region, Gene,
      SNP = !!sym(snp_col),        # 动态引用SNP列（topSNP）
      Beta = !!sym(beta_col),      # 动态引用效应量列
      SE = !!sym(se_col),          # 动态引用标准误列
      P_Value = !!sym(p_col),      # 动态引用P值列
      FDR, OR, OR_95CI,
      p_HEIDI = !!sym(p_heidi_col),# 核心修改：动态引用p_HEIDI列（解决“object not found”）
      EQTL_Prefix
    )
  
  # 5. 保存结果（不变）
  write.xlsx(
    processed_data,
    file.path(result_dir, paste0(outcome_name, "-", region_name, "_processed.xlsx")),
    row.names = FALSE
  )
  
  return(processed_data)
}

# 5. 批量运行SMR（结局×脑区，双重循环）
# 存储所有结局的结果
all_smr_results <- list()

# 循环1：遍历两个结局（ADHD+癫痫）
for (outcome in outcome_list) {
  outcome_name <- outcome$name
  outcome_path <- outcome$path
  cat(paste0("\n===== 开始处理结局：", outcome_name, " ====="))
  
  # 循环2：遍历13个脑区
  region_results <- map2_dfr(
    brain_regions$region_name,
    brain_regions$eqtl_prefix,
    function(region, eqtl) {
      run_smr(
        outcome_name = outcome_name,
        outcome_path = outcome_path,
        region_name = region,
        eqtl_prefix = eqtl
      )
    }
  )
  
  # 记录当前结局的结果
  if (!is.null(region_results) && nrow(region_results) > 0) {
    all_smr_results[[outcome_name]] <- region_results
    cat(paste0("\n✅ 结局", outcome_name, "处理完成，共", nrow(region_results), "条记录\n"))
  } else {
    all_smr_results[[outcome_name]] <- data.frame()
    cat(paste0("\n⚠️  结局", outcome_name, "无有效结果\n"))
  }
}

# 6. 合并所有结局的结果
if (length(all_smr_results) == 0) stop("❌ 所有结局均无有效结果，终止分析")

# 合并ADHD和癫痫的结果（按Gene+Brain_Region匹配）
merged_results <- all_smr_results[["ADHD"]] %>%
  rename_with(~paste0(., "_ADHD"), c("Beta", "SE", "P_Value", "FDR", "OR", "OR_95CI", "p_HEIDI")) %>%
  left_join(
    all_smr_results[["Epilepsy"]] %>%
      rename_with(~paste0(., "_Epilepsy"), c("Beta", "SE", "P_Value", "FDR", "OR", "OR_95CI", "p_HEIDI")) %>%
      dplyr::select(Gene, Brain_Region, ends_with("_Epilepsy")),
    by = c("Gene", "Brain_Region")
  ) %>%
  filter(!is.na(P_Value_Epilepsy))

cat(paste0("\n✅ 合并后共", nrow(merged_results), "个ADHD-癫痫共享基因-脑区组合\n"))

# 7. 筛选显著共享基因（与文章Figure 6标准一致）
significant_genes <- merged_results %>%
  dplyr::filter(
    # 筛选条件：FDR<0.05 + HEIDI p>0.05（排除水平多效性）
    FDR_ADHD < 0.05,
    p_HEIDI_ADHD > 0.05,
    FDR_Epilepsy < 0.05,
    p_HEIDI_Epilepsy > 0.05
  ) %>%
  # 简化脑区名称（避免绘图重叠）
  dplyr::mutate(
    Region_Simplified = case_when(
      str_detect(Brain_Region, "Amygdala") ~ "Amygdala",
      str_detect(Brain_Region, "Anterior_cingulate_cortex_BA24") ~ "Ant. Cingulate (BA24)",
      str_detect(Brain_Region, "Caudate_basal_ganglia") ~ "Caudate (BG)",
      str_detect(Brain_Region, "Cerebellar_Hemisphere") ~ "Cerebellar Hem",
      str_detect(Brain_Region, "Cerebellum") ~ "Cerebellum",
      str_detect(Brain_Region, "Cortex") ~ "Cortex",
      str_detect(Brain_Region, "Frontal_Cortex_BA9") ~ "Frontal Cortex (BA9)",
      str_detect(Brain_Region, "Hippocampus") ~ "Hippocampus",
      str_detect(Brain_Region, "Hypothalamus") ~ "Hypothalamus",
      str_detect(Brain_Region, "Nucleus_accumbens_basal_ganglia") ~ "Nucleus Accumbens (BG)",
      str_detect(Brain_Region, "Putamen_basal_ganglia") ~ "Putamen (BG)",
      str_detect(Brain_Region, "Spinal_cord_cervical_c-1") ~ "Spinal Cord (C1)",
      str_detect(Brain_Region, "Substantia_nigra") ~ "Substantia Nigra",
      TRUE ~ Brain_Region
    ),
    # 计算平均显著性（用于绘图颜色）
    Mean_NegLog10P = (-log10(P_Value_ADHD) - log10(P_Value_Epilepsy))/2
  ) %>%
  distinct(Gene, Region_Simplified, .keep_all = TRUE)

# 容错：无显著基因时的处理
if (nrow(significant_genes) == 0) {
  warning("⚠️ 无符合条件的显著共享基因（FDR<0.05 + HEIDI>0.05）")
  # 创建空数据框避免后续绘图崩溃
  significant_genes <- data.frame(
    Gene = character(),
    Region_Simplified = character(),
    Mean_NegLog10P = numeric(),
    FDR_ADHD = numeric(),
    FDR_Epilepsy = numeric(),
    OR_ADHD = numeric(),
    OR_Epilepsy = numeric(),
    OR_95CI_ADHD = character(),
    OR_95CI_Epilepsy = character(),
    stringsAsFactors = FALSE
  )
} else {
  cat(paste0("\n✅ 筛选到", nrow(significant_genes), "个显著共享基因-脑区组合（共享基因数：", length(unique(significant_genes$Gene)), "）\n"))
  # 保存显著基因结果（示例代码：Excel格式）
  write.xlsx(
    significant_genes %>%
      select(
        Gene,
        Brain_Region = Region_Simplified,
        Mean_Significance = Mean_NegLog10P,
        FDR_ADHD, FDR_Epilepsy,
        OR_ADHD = OR_ADHD, OR_95CI_ADHD,
        OR_Epilepsy = OR_Epilepsy, OR_95CI_Epilepsy,
        p_HEIDI_ADHD, p_HEIDI_Epilepsy
      ),
    file.path(result_dir, "ADHD_Epilepsy_Significant_Shared_Genes.xlsx"),
    row.names = FALSE
  )
}

# 8. 绘制Figure 6样式图（脑区×基因关联图）
if (nrow(significant_genes) > 0) {
  figure6_plot <- ggplot(
    significant_genes,
    aes(
      x = Region_Simplified,
      y = Gene,
      color = Mean_NegLog10P,
      size = Mean_NegLog10P
    )
  ) +
    geom_point(shape = 16, alpha = 0.8) +
    # 颜色映射（青色系，贴近文章风格）
    scale_color_gradient(
      low = "#4ECDC4",
      high = "#1A535C",
      name = "Mean -log10(p)",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    # 点大小（显著性越高越大）
    scale_size_range(range = c(3, 6), name = "Mean -log10(p)") +
    # 主题（示例代码：简洁无网格）
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),  # 无网格线（示例代码风格）
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    # 标签（匹配数据来源）
    labs(
      x = "GTEx V8 Brain Regions (13个脑区)",
      y = "Shared Functional Genes (ADHD-Epilepsy)",
      title = "Shared Genes of ADHD and Epilepsy in GTEx Brain Regions (SMR Analysis)"
    ) +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))
  
  # 保存图片（存子文件夹）
  ggsave(
    file.path(result_dir, "ADHD_Epilepsy_Figure6_GTEx13.pdf"),
    plot = figure6_plot,
    width = 14,
    height = 8 + length(unique(significant_genes$Gene)) * 0.2,
    dpi = 300,
    device = "pdf"
  )
  cat(paste0("\n📊 Figure 6已保存至：", file.path(result_dir, "ADHD_Epilepsy_Figure6_GTEx13.pdf"), "\n"))
}

# 9. 分析总结（示例代码：统计报告）
summary_report <- data.frame(
  Outcome = sapply(all_smr_results, function(x) if (nrow(x) > 0) unique(x$Outcome) else NA),
  Total_Records = sapply(all_smr_results, nrow),
  Total_Brain_Regions = sapply(all_smr_results, function(x) if (nrow(x) > 0) length(unique(x$Brain_Region)) else 0),
  Significant_Shared_Genes = if (nrow(significant_genes) > 0) length(unique(significant_genes$Gene)) else 0,
  Significant_Combinations = nrow(significant_genes),
  Result_Directory = result_dir,
  stringsAsFactors = FALSE
) %>% filter(!is.na(Outcome))

# 保存总结报告（示例代码：Excel统计）
write.xlsx(
  summary_report,
  file.path(workpath, "SMR_Analysis_Summary_Report.xlsx"),
  row.names = FALSE
)

# 10. 完成提示
cat("\n✅ 全部SMR分析流程完成！\n")
cat("📁 所有结果文件存于：", result_dir, "\n")
cat("📄 分析总结报告：", file.path(workpath, "SMR_Analysis_Summary_Report.xlsx"), "\n")
if (nrow(significant_genes) > 0) {
  cat("📋 显著共享基因表：", file.path(result_dir, "ADHD_Epilepsy_Significant_Shared_Genes.xlsx"), "\n")
}