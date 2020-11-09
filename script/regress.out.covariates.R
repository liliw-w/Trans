# dir.data is the directory where expression and covariates files are stored.
# Files should include:
#  1. ex.rdata: expression data stored in rdata form, including objects ex and gnames.
#  2. covariates.txt: covariates files for expression and snps, with header, no rownames, \t separated.


rm(list = ls())

#parameter = commandArgs(trailingOnly = T)
#dir.data = parameter[1] #"/project2/xuanyao/llw/breastcancerTCGA/txt_file/"
#paste0(dir.data, "covariates.txt")


extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}


# import expression and covariates data
load(snakemake@input[["file_ex"]])
colnames(ex) = gnames
cov_all = read.table(snakemake@input[["file_covariates"]],
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)

# regress out covariates
ex_cov_regressed = apply(ex, 2, function(y) extract_residual(y, as.matrix(cov_all)))
saveRDS(ex_cov_regressed, snakemake@output[["file_ex_var_regressed"]])
