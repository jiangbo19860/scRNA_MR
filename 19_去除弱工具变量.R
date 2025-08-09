# 对经过连锁不平衡（LD）聚类后的 GWAS 暴露数据进行进一步处理，补充关键遗传参数（等位基因频率）并筛选出强工具变量。输入是exposure_clumped.csv。2个输出：exposure.eaf.csv和exposure.F.csv。

rm(list = ls())  # 清空工作空间
## 安装必要的包
##install.packages("devtools")
#devtools::install_github("mrcieu/ieugwasr", force = TRUE)
pacman::p_load(
  ieugwasr,  # 用于访问IEU GWAS数据库
  httr,      # 用于HTTP请求
  jsonlite,   # 用于处理JSON数据
  here
)
here()

# 读取输入文件
inputFile = here("3_outputs/exposure_clumped.csv")
dat = read.csv(inputFile, header = TRUE, sep = ",", check.names = FALSE)

# 定义函数获取每个SNP的eaf值（等位基因频率）
# 通过访问外部数据库（Ensembl）获得每个SNP的效应等位基因频率eaf值，函数名为SNP2eaf，参数为数据框dat，目的是遍历数据框中的每一行，获取每个SNP的信息
SNP2eaf = function(dat) {
  for (i in 1:nrow(dat)) {
    # 输出进度信息
    if (i %% 10 == 0) {
      print(paste0(i, " SNP is finished!"))
    }
    # 获取当前SNP及其等位基因信息
    snpID = dat[i, "SNP"]
    localEffectAllele = dat[i, "effect_allele.exposure"]
    localOtherAllele = dat[i, "other_allele.exposure"]
    # 构建访问Ensembl数据库的URL
    website = paste0("http://rest.ensembl.org/variation/Homo_sapiens/", snpID, "?content-type=application/json;pops=1")
    # 通过API获取SNP信息
    info = httr::content(httr::GET(website))
    snpInfo = jsonlite::fromJSON(jsonlite::toJSON(info))$populations
    # 筛选出1000 Genomes Phase 3项目中欧洲人群的频率信息
    snpInfo = snpInfo[snpInfo$population == "1000GENOMES:phase_3:EUR", ]
    webEffectAllele = ifelse(is.null(snpInfo[1, "allele"][[1]]), "unknown", snpInfo[1, "allele"][[1]])
    webOtherAllele = ifelse(is.null(snpInfo[2, "allele"][[1]]), "unknown", snpInfo[2, "allele"][[1]])
    webEaf = snpInfo[1, "frequency"][[1]]
    # 比较本地和网络获取的等位基因信息，判断是否一致，设置eaf值
    if ((webEffectAllele == localEffectAllele) & (webOtherAllele == localOtherAllele)) {
      dat$eaf.exposure[i] = webEaf
    } else if ((webEffectAllele == localOtherAllele) & (webOtherAllele == localEffectAllele)) {
      dat$eaf.exposure[i] = 1 - webEaf
    } else {
      revWebEffectAllele = chartr("CGAT", "GCTA", webEffectAllele)
      revWebOtherAllele = chartr("CGAT", "GCTA", webOtherAllele)
      if ((revWebEffectAllele == localEffectAllele) & (revWebOtherAllele == localOtherAllele)) {
        dat$eaf.exposure[i] = webEaf
      } else if ((revWebEffectAllele == localOtherAllele) & (revWebOtherAllele == localEffectAllele)) {
        dat$eaf.exposure[i] = 1 - webEaf
      } else {
        dat$eaf.exposure[i] = ifelse(localEffectAllele == webEffectAllele, webEaf, 1 - webEaf)
      }
    }
  }
  return(dat)
}

# 调用函数获取eaf值
dat = SNP2eaf(dat)
# 将结果写入文件
write.csv(dat, file = "3_outputs/exposure.eaf.csv", row.names = FALSE)

# 计算F值 --------------------------------------------------------------------
# R²（决定系数）是SNP 能够解释暴露表型变异的比例，取值范围 0~1，值越大说明 SNP 与暴露的关联越强。公式的分子和分母分别代表 “遗传效应带来的变异”（遗传效应方差） 和 “总变异”（总表型方差）：beta：SNP 对暴露的遗传效应值（从 GWAS 中获得）；eaf：效应等位基因频率。
dat$R2 <- (2 * dat$beta.exposure * dat$beta.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure) /
             (2 * dat$beta.exposure * dat$beta.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure) +
                2 * dat$se.exposure * dat$se.exposure * dat$samplesize.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure)))

# F值：评估工具变量的强度，F > 10 通常被认为是强工具变量。公式：F = R² * (n - k - 1) / (1 - R²)，其中 n 是样本量，k 是工具变量的数量（这里为 2，因为每个 SNP 有效等位基因和其他等位基因）。(N - 2)：自由度调整（N 为样本量，减 2 是因为回归模型中估计了 2 个参数：截距和效应值 beta）；F 值与 R² 正相关（R² 越大，F 值越大），与样本量正相关（样本量越大，相同 R² 下 F 值越大）。
dat$F <- dat$R2 * (dat$samplesize.exposure - 2) / (1 - dat$R2)

# 过滤F值大于10的数据，并写入文件file = "3_outputs/exposure.eaf.csv", row.names = FALSE)
outTab = dat[as.numeric(dat$F) > 10, ]
write.csv(outTab, file = "3_outputs/exposure.F.csv", row.names = FALSE)
