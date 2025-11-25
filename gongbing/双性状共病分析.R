
library(QTLMR)

# PMID: 38425786（流程一致，结果无法完全一致）

#########################数据转换###################################

trait1 <- "Frailty"                      #性状1的名称，需要指定该参数

trait2 <- "Insomnia"                     #性状2的名称，需要指定该参数

trait1_samplesize <- 175226              #性状1的总样本量，需要指定该参数

trait2_samplesize <- 462341              #性状2的总样本量，需要指定该参数


format_data_Catalog(GWASfile = "./1.data/GCST90020053_buildGRCh37.tsv",
                    GWAS_name = trait1,
                    type = "exposure",
                    samplesize = trait1_samplesize,
                    get_SNP = FALSE,
                    build = NULL,
                    save_path = "./2.pre_data",
                    Twosample_dat = TRUE,
                    SMR_dat = TRUE,
                    MTAG_dat = TRUE,
                    METAL_dat = TRUE)


format_data_IEU_VCF(IEU_VCF = "./1.data/ukb-b-3957.vcf.gz",
                    type = "outcome",
                    min_pval = 0,
                    samplesize = trait2_samplesize,
                    Twosample_dat = TRUE,
                    SMR_dat = TRUE,
                    MTAG_dat = TRUE,
                    METAL_dat = TRUE,
                    save_name = trait2,
                    save_path = "./2.pre_data")

#根据上面函数生成的格式数据，自行指定修改性状1与性状2的TwosampleMR与MTAG格式数据名称

TwosampleMR_dat_trait1 <- "Frailty_TwosampleMR.txt"

TwosampleMR_dat_trait2 <- "Insomnia_TwosampleMR.txt"

MTAG_dat_trait1 <- "Frailty_MTAG.txt"

MTAG_dat_trait2 <- "Insomnia_MTAG.txt"



#########################LDSC分析###################################

LDSC_py_sumstats(GWASfile = paste0("./2.pre_data/",TwosampleMR_dat_trait1),
                 GWAS_name = trait1,
                 N = NULL,
                 merge_alleles = "./LDSC/w_hm3.noMHC.snplist",
                 save_path = paste0("./LDSC/",trait1))

LDSC_py_sumstats(GWASfile = paste0("./2.pre_data/",TwosampleMR_dat_trait2),
                 GWAS_name = trait2,
                 N = NULL,
                 merge_alleles = "./LDSC/w_hm3.noMHC.snplist",
                 save_path = paste0("./LDSC/",trait2))


#考虑截距存在
LDSC_py_rg(test_help = FALSE,
           Sumstatsfile = c(paste0("./LDSC/",trait1,"/",trait1,".sumstats.gz"),
                            paste0("./LDSC/",trait2,"/",trait2,".sumstats.gz")),
           ref_ld_chr = "./LDSC/eur_w_ld_chr",
           w_ld_chr = "./LDSC/eur_w_ld_chr",
           opt_arguments = NULL,
           save_name = "rg",
           save_path = "./LDSC/LDSC_rg")

#不考虑截距存在，在opt_arguments中，将--no-intercept参数写入即可。
LDSC_py_rg(test_help = FALSE,
           Sumstatsfile = c(paste0("./LDSC/",trait1,"/",trait1,".sumstats.gz"),
                            paste0("./LDSC/",trait2,"/",trait2,".sumstats.gz")),
           ref_ld_chr = "./LDSC/eur_w_ld_chr",
           w_ld_chr = "./LDSC/eur_w_ld_chr",
           opt_arguments = "--no-intercept",
           save_name = "no_rg",
           save_path = "./LDSC/LDSC_rg")


#########################S_LDSC分析#################################

#必须将工作路径设置在与"GTEx.ldcts"文件在相同路径下!

setwd("./LDSC")

LDSC_py_ph(test_help = FALSE,
           Sumstatsfile = paste0("./",trait1,"/",trait1,".sumstats.gz"),
           ref_ld_chr = "./LDSCORE_1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.",
           w_ld_chr = "./weights_hm3_no_hla/weights.",
           frqfile_chr = "./1000G_Phase3_frq/1000G.EUR.QC.",
           opt_arguments = NULL,
           save_name = trait1,
           save_path = "./S_LDSC")

LDSC_py_ph(test_help = FALSE,
           Sumstatsfile = paste0("./",trait2,"/",trait2,".sumstats.gz"),
           ref_ld_chr = "./LDSCORE_1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.",
           w_ld_chr = "./weights_hm3_no_hla/weights.",
           frqfile_chr = "./1000G_Phase3_frq/1000G.EUR.QC.",
           opt_arguments = NULL,
           save_name = trait2,
           save_path = "./S_LDSC")


S_LDSC_plot <- Visualizing_bidir_plot(trait1_file = paste0("./S_LDSC/",trait1,".results"),
                                      trait2_file = paste0("./S_LDSC/",trait2,".results"),
                                      Pval_threshold = 0.05,
                                      Pval_column = "Enrichment_p",
                                      category_column = "Category",
                                      trait1_name = trait1,
                                      trait2_name = trait2,
                                      color_trait1 = "#5A9BD4",
                                      color_trait2 = "#E6958A",
                                      y_axis_label = "Category",
                                      y_axis_label_size = 10,
                                      y_axis_font_size = 8,
                                      x_axis_font_size = 8,
                                      x_axis_label = expression(-log[10](Enrichment_pval)),
                                      x_axis_label_size = 10,
                                      x_axis_breaks = c(-16, -12, -8, -4, 0, 4, 8, 12, 16),
                                      legend_title_size = 12,
                                      legend_text_size = 10,
                                      pdf_width = 10,
                                      pdf_height = 5.5,
                                      save_name = "s_ldsc",
                                      save_path = "./S_LDSC/")

########################LDSC组织特异性分析#########################

LDSC_py_SEG(test_help = FALSE,
            Sumstatsfile = paste0("./",trait1,"/",trait1,".sumstats.gz"),
            ref_ld_chr = "./1000G_EUR_Phase3_baseline/baseline.",
            w_ld_chr = "./weights_hm3_no_hla/weights.",
            ref_ld_chr_cts = "GTEx.ldcts",
            opt_arguments = NULL,
            save_name = trait1,
            save_path = "./LDSC_SEG")

LDSC_py_SEG(test_help = FALSE,
            Sumstatsfile = paste0("./",trait2,"/",trait2,".sumstats.gz"),
            ref_ld_chr = "./1000G_EUR_Phase3_baseline/baseline.",
            w_ld_chr = "./weights_hm3_no_hla/weights.",
            ref_ld_chr_cts = "GTEx.ldcts",
            opt_arguments = NULL,
            save_name = trait2,
            save_path = "./LDSC_SEG")

SEG_plot <- Visualizing_bidir_plot(trait1_file = paste0("./LDSC_SEG/",trait1,".cell_type_results.txt"),
                                   trait2_file = paste0("./LDSC_SEG/",trait2,".cell_type_results.txt"),
                                   Pval_threshold = 1,
                                   Pval_column = "Coefficient_P_value",
                                   category_column = "Name",
                                   trait1_name = trait1,
                                   trait2_name = trait2,
                                   color_trait1 = "#5A9BD4",
                                   color_trait2 = "#E6958A",
                                   y_axis_label = "Tissue",
                                   y_axis_label_size = 10,
                                   y_axis_font_size = 8,
                                   x_axis_font_size = 8,
                                   x_axis_label = expression(-log[10](Coefficient_P_value)),
                                   x_axis_label_size = 10,
                                   x_axis_breaks = c(-6, -4, -2, 0, 2, 4, 6),
                                   legend_title_size = 12,
                                   legend_text_size = 10,
                                   pdf_width = 10,
                                   pdf_height = 7,
                                   save_name = "LDSC_SEG",
                                   save_path = "./LDSC_SEG")


#########################MAGMA分析###################################

setwd("../") #返回至上一级作为工作目录


# 对性状1进行组织特异性分析



gene_based_trait1 <- MAGMA_gene_based(GWAS_file = paste0("./2.pre_data/",MTAG_dat_trait1),
                                      bfile_1000G = "./MAGMA/g1000_eur/g1000_eur",
                                      gene_loc = "./MAGMA/ENSGv110.coding.genes.txt",
                                      set_annot = "./MAGMA/MSigDB_20231Hs_MAGMA.txt",
                                      SNP_P_col = c(3, 9),
                                      samplesize_col = "N",
                                      save_name = trait1,
                                      save_path = paste0("./MAGMA/",trait1))

Tissue_specific_trait1 <- MAGMA_Tissue_specific(genes_raw = paste0("./MAGMA/",trait1,"/",trait1,".genes.raw"),
                                                gene_covar = "./MAGMA/gtex_v8_ts_avg_log2TPM.txt",
                                                save_name = trait1,
                                                save_path = paste0("./MAGMA/",trait1))

write.csv(Tissue_specific_trait1,file = paste0("./MAGMA/",trait1,"/Tissue_specific_",trait1,".csv"),row.names = F)


# 对性状2进行组织特异性分析

gene_based_trait2 <- MAGMA_gene_based(GWAS_file = paste0("./2.pre_data/",MTAG_dat_trait2),
                                      bfile_1000G = "./MAGMA/g1000_eur/g1000_eur",
                                      gene_loc = "./MAGMA/ENSGv110.coding.genes.txt",
                                      set_annot = "./MAGMA/MSigDB_20231Hs_MAGMA.txt",
                                      SNP_P_col = c(3, 9),
                                      samplesize_col = "N",
                                      save_name = trait2,
                                      save_path = paste0("./MAGMA/",trait2))

Tissue_specific_trait2 <- MAGMA_Tissue_specific(genes_raw = paste0("./MAGMA/",trait2,"/",trait2,".genes.raw"),
                                                gene_covar = "./MAGMA/gtex_v8_ts_avg_log2TPM.txt",
                                                save_name = trait2,
                                                save_path = paste0("./MAGMA/",trait2))

write.csv(Tissue_specific_trait2,file = paste0("./MAGMA/",trait2,"/Tissue_specific_",trait2,".csv"),row.names = F)


SEG_plot <- Visualizing_bidir_plot(trait1_file = paste0("./MAGMA/",trait1,"/Tissue_specific_",trait1,".csv"),
                                   trait2_file = paste0("./MAGMA/",trait2,"/Tissue_specific_",trait2,".csv"),
                                   Pval_threshold = 1,
                                   Pval_column = "P",
                                   category_column = "FULL_NAME",
                                   trait1_name = trait1,
                                   trait2_name = trait2,
                                   color_trait1 = "#5A9BD4",
                                   color_trait2 = "#E6958A",
                                   y_axis_label = "Tissue",
                                   y_axis_label_size = 10,
                                   y_axis_font_size = 8,
                                   x_axis_font_size = 8,
                                   x_axis_label = expression(-log[10](P)),
                                   x_axis_label_size = 10,
                                   x_axis_breaks = c(-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12),
                                   legend_title_size = 12,
                                   legend_text_size = 10,
                                   pdf_width = 10,
                                   pdf_height = 7,
                                   save_name = "MAGMA",
                                   save_path = "./MAGMA")


#########################HDL分析###################################

#
# dat <- HDL_rg_parallel(GWASfile = c(paste0("./2.pre_data/",MTAG_dat_trait1),
#                                     paste0("./2.pre_data/",MTAG_dat_trait2)),
#                        LD.path = "G:/HDL/UKB_imputed_SVD_eigen99_extraction",   #数据下载地址 https://github.com/zhenin/HDL/wiki/FAQ
#                        Nref = 335265,
#                        N0 = NULL,
#                        numCores = 10,
#                        eigen.cut = "automatic",
#                        jackknife.df = FALSE,
#                        intercept.output = FALSE,
#                        fill.missing.N = NULL,
#                        lim = exp(-18),
#                        verbose = FALSE,
#                        save_name = "rg",
#                        save_path = "./HDL")

#########################HESS分析###################################



HESS_format_data(GWASfile = paste0("./2.pre_data/",TwosampleMR_dat_trait1),
                 N = NULL,
                 FOCUS_munge = FALSE,
                 gwas_loci = TRUE,
                 save_name = trait1,
                 save_path = paste0("./HESS/",trait1))

HESS_format_data(GWASfile = paste0("./2.pre_data/",TwosampleMR_dat_trait2),
                 N = NULL,
                 FOCUS_munge = FALSE,
                 gwas_loci = TRUE,
                 save_name = trait2,
                 save_path = paste0("./HESS/",trait2))


dat_HESS <- HESS_estimate_correlation(Sumstatsfile = c(paste0("./HESS/",trait1,"/",trait1,".sumstats"),
                                                       paste0("./HESS/",trait2,"/",trait2,".sumstats")),
                                      bfile = "./HESS/1kg_eur_1pct/1kg_eur_1pct_chr",
                                      partition = "./HESS/ldetect-data/EUR/fourier_ls-all.bed",
                                      reinflate_lambda_gc = NULL,
                                      tot_hsqg = NULL,
                                      num_shared = 0,
                                      pheno_cor = 0,
                                      start_chr = 1,
                                      end_chr = 22,
                                      opt_arguments_step1 = NULL,
                                      opt_arguments_step2 = NULL,
                                      opt_arguments_step3 = NULL,
                                      save_name = "step3",
                                      save_path = "./HESS/HESS",
                                      cores = 1)

HESS_Visualizing_estimates(local_rhog_est = "./HESS/HESS/step3.txt",
                           local_hsqg_est = c("./HESS/HESS/step2_trait1.txt","./HESS/HESS/step2_trait2.txt"),
                           trait_names = c(trait1, trait2),
                           save_name = "Visualizing",
                           save_path = "./HESS")


HESS_res <-data.table::fread("./HESS/HESS/step3.txt",showProgress = F) %>%
  dplyr::arrange(p)

HESS_res$fdr <- p.adjust(HESS_res$p,method = "bonferroni")

HESS_res <- subset(HESS_res,HESS_res$fdr<0.05)


#########################MTAG分析###################################

dat_MTAG <- MTAG_cross_gwas(GWASfile = c(paste0("./2.pre_data/",MTAG_dat_trait1),
                                         paste0("./2.pre_data/",MTAG_dat_trait2)),
                            snp_col = "SNP",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            beta_col = "beta",
                            se_col = "se",
                            pval_col = "pval",
                            eaf_col = "eaf",
                            samplesize_col = "N",
                            Zscores_col = "Z",
                            maf_min = 0.01,
                            n_min = NULL,
                            no_overlap = FALSE,
                            perfect_gencov = FALSE,
                            equal_h2 = FALSE,
                            fdr = TRUE,
                            opt_arguments = NULL,
                            save_name = "MTAG",
                            save_path = "./MTAG")

#########################CPASSOC分析###################################

dat_CPASSOC <- CPASSOC_meta_gwas(GWASfile = c(paste0("./2.pre_data/",MTAG_dat_trait1),
                                              paste0("./2.pre_data/",MTAG_dat_trait2)),
                                 GWAS_samplesizes = c(trait1_samplesize,trait2_samplesize),
                                 GWAS_name = c(trait1, trait2),
                                 mvrnorm_n = 1e+06,
                                 methods_sig = "SHet",
                                 sig_pval = 5e-08,
                                 save_name = paste0(trait1,"_",trait2),
                                 save_path = "./CPASSOC",
                                 cores = 4)


#########################novel_SNPs分析###################################

dat_CPASSOC <- data.table::fread("./CPASSOC/cpassoc_sig.txt")[,c(1,8)]

dat_CPASSOC <- subset(dat_CPASSOC,dat_CPASSOC$p_SHet<5e-8)

dat_MTAG <- data.table::fread("./MTAG/merged_sig_pval.txt")

data <- merge(dat_MTAG,dat_CPASSOC,by = "SNP")

sig_ind_snp <- clump_data_local_Online(data,
                                       snp_col = "SNP",
                                       pval_col = "p_SHet",
                                       clump_pval = 5e-08,
                                       clump_kb = 500,
                                       clump_r2 = 0.3,
                                       bfile_1000G = "./1.data/1000G/EUR")

sig_pre_SNP <- subset(sig_ind_snp,sig_ind_snp$trait1_pval<5e-8&sig_ind_snp$trait2_pval<5e-8)

r2_SNP <- clump_LD_R2(sig_pre_SNP,
                      snp_col = "SNP",
                      clump_kb = 1000,
                      clump_r2 = 0.3,
                      bfile_1000G = "./1.data/1000G/EUR")

dat3 <- subset(sig_ind_snp ,!sig_ind_snp$SNP%in%r2_SNP)


dat4 <- dat3 %>%subset(.,trait1_pval<5e-8 & trait2_pval<5e-8)

novel_SNPs <- subset(dat3,!dat3$SNP%in%dat4$SNP)


#合并hess结果，此步可做，可不做

hess_dat <- data.table::fread("./HESS/HESS/step3.txt")
hess_dat$fdr <- p.adjust(hess_dat$p,method = "bonferroni")
hess_dat <- subset(hess_dat,hess_dat$fdr<0.05)

novel_SNPs_merge <- c()

for (j in 1:nrow(hess_dat)) {

  hess_data <- hess_dat[j,]

  novel_SNPs_new <- novel_SNPs %>% dplyr::filter((trait1_CHR == hess_data$chr)&
                                                   (trait1_BP > hess_data$start) &
                                                   (trait1_BP < hess_data$end))

  novel_SNPs_merge <- rbind(novel_SNPs_merge,novel_SNPs_new)
}


#########################MR分析#########################################

#双向GSMR分析

dat_GSMR <- GSMR2_QTLMR(exposurefile =paste0("./2.pre_data/",TwosampleMR_dat_trait1),
                        exposure_name = trait1,
                        outcomefile = paste0("./2.pre_data/",TwosampleMR_dat_trait2),
                        outcome_name = trait2,
                        bfile_1000G = "./1.data/1000G/EUR",
                        n_ref = 503,
                        gwas_thresh = 5e-08,
                        heidi_outlier_flag = TRUE,
                        single_snp_heidi_thresh = 0.01,
                        multi_snps_heidi_thresh = 0.01,
                        nsnps_thresh = 10,
                        ld_r2_thresh = 0.05,
                        ld_fdr_thresh = 0.05,
                        gsmr2_beta = TRUE,
                        GSMR_plot = TRUE,
                        width = 6,
                        height = 6,
                        suppressWarnings = TRUE,
                        Bi_GSMR = TRUE,
                        save_name = "GSMR2",
                        save_path = "./GSMR2",
                        thread_num = 2)

write.csv(dat_GSMR,file = "./GSMR2/dat_GSMR.csv",,row.names = F)


#正向MR分析

library(TwoSampleMR)

exposure_dat <- data.table::fread(paste0("./2.pre_data/",TwosampleMR_dat_trait1))

exposure_dat <- clump_data_local_Online(exposure_dat,
                                        snp_col = "SNP",
                                        pval_col = "pval.exposure",
                                        clump_method = "PVAL",
                                        beta_col = "beta.exposure",
                                        se_col = "se.exposure",
                                        clump_pval = 5e-08,
                                        clump_kb = 10000,
                                        clump_r2 = 0.001,
                                        pop = "EUR",
                                        bfile_1000G = "./1.data/1000G/EUR")

outcome_dat <- as.data.frame(data.table::fread(paste0("./2.pre_data/",TwosampleMR_dat_trait2)))

dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)

res1 <- xQTL_mr(dat, FDR_method = NULL, PVE = F)

res1$exposure <- trait1

res1$outcome <- trait2

#反向MR分析

exposure_dat <- as.data.frame(data.table::fread(paste0("./2.pre_data/",TwosampleMR_dat_trait2)))

if (ncol(exposure_dat[, grep(pattern = "outcome", x = colnames(exposure_dat))]) > 0) {
  colnames(exposure_dat) = gsub(pattern = "outcome",
                                replacement = "exposure",
                                x = colnames(exposure_dat))}

exposure_dat <- clump_data_local_Online(exposure_dat,
                                        snp_col = "SNP",
                                        pval_col = "pval.exposure",
                                        clump_method = "PVAL",
                                        beta_col = "beta.exposure",
                                        se_col = "se.exposure",
                                        clump_pval = 5e-08,
                                        clump_kb = 10000,
                                        clump_r2 = 0.001,
                                        pop = "EUR",
                                        bfile_1000G = "./1.data/1000G/EUR")


outcome_dat <- as.data.frame(data.table::fread(paste0("./2.pre_data/",TwosampleMR_dat_trait1)))

if (ncol(outcome_dat[, grep(pattern = "exposure", x = colnames(outcome_dat))]) > 0) {
  colnames(outcome_dat) = gsub(pattern = "exposure",
                               replacement = "outcome",
                               x = colnames(outcome_dat))}

dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)

res2 <- xQTL_mr(dat, FDR_method = NULL, PVE = F)

res2$exposure <- trait2

res2$outcome <- trait1


res <- rbind(res1,res2,dat_GSMR)%>%
  dplyr::mutate(method = dplyr::recode(method, "Inverse variance weighted" = "IVW"))%>%
  dplyr::arrange(desc(exposure))

dir.create("./MR")

save(res,file = "./MR/res.Rdata")

Visualizing_MR_forest(res,
                      fdr_col = "fdr",
                      ci_col = "#7A378B",
                      ci_fill = "#008B00",
                      ci_lwd = 0.6,
                      ci_sizes = 0.3,
                      xlim = c(0,3),
                      ticks_at = c(0,1,2,3),
                      space_padding = 20,
                      arrange_OR = T,
                      save_plot = T,
                      plot_pdf = "forest_plot",
                      width = 8,
                      height = 5,
                      save_path = "./MR")


#########################SMR分析###################################

data <- as.data.frame(data.table::fread("./2.pre_data/Frailty_SMR_Catalog.ma"))

data <- subset(data,data$p<5e-8)

data.table::fwrite(x = as.data.frame(data), file = "./2.pre_data/Frailty_SMR_Catalog_sig.ma",
                   row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


data <- as.data.frame(data.table::fread("./2.pre_data/ukb-b-3957_SMR_IEU.ma"))

data <- subset(data,data$p<5e-8)

data.table::fwrite(x = as.data.frame(data), file = "./2.pre_data/ukb-b-3957_SMR_IEU_sig.ma",
                   row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# 分别将Frailty_SMR_Catalog_sig.ma与ukb-b-3957_SMR_IEU_sig.ma上传到SMR在线网站，对共同显著的组织进行SMR分析
# https://yanglab.westlake.edu.cn/smr-portal/
