rm(list = ls())

if(!require("TwoSampleMR")){
  library(devtools)
  install_github("MRCIEU/TwoSampleMR")
}


library(TwoSampleMR)

pop = "EUR"
gwasPhenocode = "30090"


file_qtlMaf = "/project2/xuanyao/llw/DGN/data/chr_all.frq"

file_qtlColocReg = "/scratch/midway2/liliw1/coloc/data/qtlColocReg.txt.gz"
file_gwas_bgz = list.files(path = "/scratch/midway2/liliw1/UKB_nealelab", pattern = paste0("^continuous.", gwasPhenocode, ".*.tsv.bgz$"), full.names = TRUE)
file_snp_meta_qtl = "/scratch/midway2/liliw1/coloc/snp.meta.qtl.txt.gz"


qtlColocReg = fread(file_qtlColocReg)
snp_meta_qtl = fread(file_snp_meta_qtl)

gwasCol = c("chr", "pos", "ref", "alt", paste(c("af", "beta", "se", "pval", "low_confidence"), pop, sep = "_") )
gwas = fread(cmd = paste("gunzip -c", file_gwas_bgz), select = gwasCol)
gwas$if_good = !is.na(gwas[[paste0("pval_", pop)]]) & !gwas[[paste0("low_confidence_", pop)]]

gwas = gwas[gwas$if_good, ]
gwas$low_confidence_EUR = NULL
gwas$if_good = NULL

a = paste(gwas$chr, gwas$pos, sep = ":")
str(a)
gwas$SNP_ID = a

gwas$SNP_ID = paste(gwas$chr, gwas$pos, sep = ":")
gwas = left_join(x = gwas, y = snp_meta_qtl %>% select(c(SNP_ID, rsid)), by = c("SNP_ID"))

gwas = gwas[!is.na(gwas$rsid), ]

fwrite(gwas, "gwas.mr.txt.gz", quote=FALSE, sep='\t')




library(TwoSampleMR)

gwasPhenocode = "30090"
tmp_file = list.files("~/Downloads", "^pheno30090.*.reg.coloc.png$")


qtlColocReg = fread("~/Downloads/qtlColocReg.txt.gz")
qtlMaf = fread("~/Downloads/chr_all.frq")
gwas = fread("~/Downloads/gwas.mr.txt.gz")

qtlMaf = qtlMaf[!duplicated(qtlMaf$SNP), ]
#if( !all(qtlColocReg$SNP_ID %in% qtlMaf$SNP) ) stop("Not all SNPs have maf!")



for(tmp in tmp_file){
  tmp_split = strsplit(tmp, ".", fixed = TRUE)[[1]]
  reg = tmp_split[2]

  qtlColocReg_reg = qtlColocReg %>% filter(Region == reg)

  qtlColocReg_reg = qtlColocReg_reg %>%
    left_join(y = qtlMaf %>% mutate(MAF = NULL), by = c("SNP_ID" = "SNP", "Chr" = "CHR")) %>%
    mutate( NCHROBS = NULL)
  qtlColocReg_reg = qtlColocReg_reg %>%
    mutate(z = sqrt(qchisq(1 - Pval, 1)) ) %>%
    mutate(beta = z, se = 1)
  qtlColocReg_reg = qtlColocReg_reg[qtlColocReg_reg$rsid != "", ]

  E <- format_data(
    qtlColocReg_reg,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "Module",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "MAF",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    pval_col = "Pval",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    #samplesize_col = "samplesize",
    #gene_col = "gene",
    #id_col = "id",
    #min_pval = 1e-200,
    z_col = "z",
    #info_col = "info",
    chr_col = "Chr",
    pos_col = "Pos",
    #log_pval = FALSE
  )
  E$samplesize = 913


  E <- clump_data(
    E,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR"
  )



  gwas$beta_EUR = abs(gwas$beta_EUR)
  O = format_data(
    gwas,
    type = 'outcome',
    snps = E$SNP,
    header = TRUE,
    #phenotype_col = "Module",
    snp_col = "rsid",
    beta_col = "beta_EUR",
    se_col = "se_EUR",
    eaf_col = "af_EUR",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    pval_col = "pval_EUR",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    #samplesize_col = "samplesize",
    #gene_col = "gene",
    #id_col = "id",
    #min_pval = 1e-200,
    #z_col = "z",
    #info_col = "info",
    chr_col = "chr",
    pos_col = "pos",
    #log_pval = FALSE
  )
  O$outcome = gwasPhenocode

  dat = harmonise_data(E, O)
  dat$beta.exposure = abs(dat$beta.exposure)
  dat$beta.outcome = abs(dat$beta.outcome)

  res = mr(dat)

  resh = mr_heterogeneity(dat)

  #list(data = dat, mr = res, heterogeneity = resh)

  p1 <- mr_scatter_plot(res, dat)
  fig <- p1[[1]]
  fig <- fig  + labs(title = paste0("pheno", gwasPhenocode, ".", reg) ) +
    scale_color_discrete(labels=sort(paste(res$method, format(res$pval, digits = 3), res$nsnp, sep = ";")) )


  ggsave(filename =  paste0("pheno", gwasPhenocode, ".", reg, ".mr.png"), fig)

  cat("Region", reg, "is done. \n")

}







