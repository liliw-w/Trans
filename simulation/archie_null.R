##############################################
############ archie null ############
##############################################
# load packages -----
rm(list = ls())
library(mvtnorm)
library(optparse)


source("/project2/xuanyao/llw/compare_to_archie/ARCHIE/codes/helper.R")
source("/project2/xuanyao/llw/compare_to_archie/ARCHIE/codes/main.R")


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    /project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds
    500
    1000000
    archie_null_cc_value_z.rds
    archie_null_selected_gene_z.rds
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_Sigma <- args[1]
N <- as.numeric(args[2])
n_sim <- as.numeric(args[3])

file_cc_null <- args[4]
file_selected_gene_null <- args[5]



# read files -----
Sigma <- as.matrix(readRDS(file_Sigma))
K <- dim(Sigma)[1]



# run archie test -----
cc_null <- array(
  dim = c(1, n_sim),
  dimnames = list(c("archie"), NULL)
)
selected_gene_null <- array(
  dim = c(K, n_sim),
  dimnames = list(NULL, NULL)
)
for(i in 1:n_sim){
  ## null z
  z_null <- rmvnorm(1, rep(0, K), Sigma) %>% t()
  
  ## archie input
  Sigma_EE <- as.matrix(1)
  Sigma_GG <- Sigma
  Sigma_GE <- z_null / sqrt(N)
  
  
  ## archie
  res.alt <- archie_work(Sigma_GE, Sigma_GG, Sigma_EE, verbose = FALSE)
  
  
  ## save
  cc_null['archie', i] <- max(res.alt$qs)
  selected_gene_null[, i] <- res.alt$us[, which.max(res.alt$qs)]
  
  
  if(i %% 1000 == 0){
    cat("Simulation: ", i, "\n")
    
    saveRDS(cc_null, file_cc_null)
    saveRDS(selected_gene_null, file_selected_gene_null)
  }
}


# print out key message or write out -----
saveRDS(cc_null, file_cc_null)
saveRDS(selected_gene_null, file_selected_gene_null)

