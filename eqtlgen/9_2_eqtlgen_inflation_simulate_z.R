###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Simulate z-scores based on the selected sigma for running PCO ###########
###############################################################
library(mvtnorm)

### I/O
module <- as.numeric(snakemake@params[['module']])
n_SNP <- 1e+4

#file_Sigma <- list.files(path = "inflation",
#                         pattern = paste0("^Sigma_DGN_module", module, "\\_K\\d*\\.rds$"),
#                         full.names = TRUE)
file_Sigma <- list.files(path = "inflation",
                         pattern = paste0("^Sigma_DGN_module", module, "\\.rds$"),
                         full.names = TRUE)
file_Z <- paste0('inflation/Z_', basename(file_Sigma))


### read data
Sigma <- as.matrix(readRDS(file_Sigma))
K <- nrow(Sigma)

### simulate z
Z.null <- rmvnorm(n_SNP, rep(0, K), Sigma)

### save Z
saveRDS(Z.null, file_Z)

