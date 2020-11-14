# dir.data is the directory where expression and covariates files are stored.
# Files should include:
#  1. ex.rdata: expression data stored in rdata form, including objects ex and gnames.
#  2. covariates.txt: covariates files for expression and snps, with header, no rownames, \t separated.


#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#dir.data = parameter[1] #"/project2/xuanyao/llw/breastcancerTCGA/txt_file/"
#paste0(dir.data, "covariates.txt")

input1 = snakemake@input[["file_ex"]]
input2 = snakemake@input[["file_covariates"]]
output1 = snakemake@output[["file_ex_var_regressed"]]

extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}


# import expression and covariates data
load(input1)
colnames(ex) = gnames
cov_all = t(as.matrix(read.table(input2,
                     sep = "\t", header = TRUE, row.names = 1,
                     stringsAsFactors = FALSE, check.names = FALSE)))

# regress out covariates
ex_cov_regressed = apply(ex, 2, function(y) extract_residual(y, cov_all))
saveRDS(ex_cov_regressed, output1)
