###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Simulate z-scores based on the selected sigma for running PCO ###########
###############################################################
module <- as.numeric(snakemake@params[['module']])

# load packages -----
library(mvtnorm)

# I/O & paras -----
n_SNP <- 1e+4
file_Sigma <- paste0("inflation/Sigma_DGN_module", module, ".rds")

## output -----
file_Z <- paste0('inflation/Z_', basename(file_Sigma))


# read files -----
Sigma <- as.matrix(readRDS(file_Sigma))
K <- nrow(Sigma)


# simulate z -----
Z.null <- rmvnorm(n_SNP, rep(0, K), Sigma)


# print out key message or write out -----
saveRDS(Z.null, file_Z)

