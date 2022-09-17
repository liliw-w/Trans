###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Estimate the sigma based on null z scores of different numbers ###########
###############################################################
module <- as.numeric(snakemake@params[['module']])


# load packages -----
library(mvtnorm)


# I/O & paras -----
thre_p_z <- 1e-4
ratios <- c(1, 5, 10, 50, 100, 150)
file_Sigma <- paste0("inflation/Sigma_DGN_module", module, ".rds")

## output -----


# read files & important paras -----
Sigma <- as.matrix(readRDS(file_Sigma))

## module size & null snps sizes of various ratios
K <- nrow(Sigma)
n_nullz_seq <- floor(K * ratios)

## z2 cutoff of p cutoff
thre_z2 <- qchisq(1-thre_p_z, 1)


# estimated sigma on null z of different sizes -----
cat(
  'For a', K, 'dimension module, its Sigma will be estimated on \n\n',
  n_nullz_seq, "NULL z's. \n\n"
)

for(n_nullz in n_nullz_seq){
  Z_null <- matrix(ncol = K)
  while (nrow(Z_null) < n_nullz) {
    tmp_Z_null <- rmvnorm(n_nullz, rep(0, K), Sigma)
    
    altSNP_ind <- apply(tmp_Z_null^2 > thre_z2, 1, sum) > 0
    
    Z_null <- rbind(Z_null, tmp_Z_null[!altSNP_ind, ])
  }
  Z_null <- Z_null[2:n_nullz, ]
  
  Sigma_nullz <- cor(Z_null)
  
  ## write out sigma of nullz -----
  file_Sigma_nullz <- paste0('inflation/Sigma_nullz', n_nullz, "_", basename(file_Sigma))
  saveRDS(Sigma_nullz, file_Sigma_nullz)
  
  cat(n_nullz, "case is done. \n")
}


###############################################################
############### for snakemake pipeline use ###############
###############################################################
file_snm <- paste0('inflation/Sigma_nullz_', basename(file_Sigma))
saveRDS(NULL, file_snm)

