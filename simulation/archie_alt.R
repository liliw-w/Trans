##############################################
############ archie ############
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
library(data.table)
library(mvtnorm)
library(optparse)


source("ARCHIE/codes/helper.R")
source("ARCHIE/codes/main.R")


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    var
    0.1
    /project2/xuanyao/llw/DGN_no_filter_on_mappability/script/
    /project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds
    0.001
    0.3
    500
    1
    1000
    archie_cc_value_addN.rds
    archie_selected_gene_addN.rds
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

change <- args[1]
oracle_thre <- as.numeric(args[2])
dir_pco <- args[3]
file_Sigma <- args[4]

var.b <- as.numeric(args[5])
caus <- as.numeric(args[6])
N <- as.numeric(args[7])

N.sample <- as.numeric(args[8])
N.sim <- as.numeric(args[9])

file_cc_alt <- args[10]
file_selected_gene_alt <- args[11]


# read files -----
# source pco test
# source(paste0(dir_pco, "ModifiedPCOMerged.R"))
# source(paste0(dir_pco, "liu.R"))
# source(paste0(dir_pco, "liumod.R"))
# source(paste0(dir_pco, "davies.R"))
# dyn.load(paste0(dir_pco, "qfc.so"))
# source(paste0(dir_pco, "ModifiedSigmaOEstimate.R"))

Sigma <- as.matrix(readRDS(file_Sigma))
K <- dim(Sigma)[1]


# prepare eigenvalues and eigenvectors as input for three methods -----
# SigmaO <- ModifiedSigmaOEstimate(Sigma)
# eigen_res <- eigen(Sigma)
# eigen_lamb <- eigen_res$values
# eigen_vec <- eigen_res$vectors

# 
# sigma_trunc_inv <- eigen_vec[, which(eigen_lamb>oracle_thre)] %*% diag(1/eigen_lamb[eigen_lamb>oracle_thre]) %*% t(eigen_vec[, which(eigen_lamb>oracle_thre)])


# p values of alt z's under various N -----
if(change == "N"){
  cat(paste0("Various ", change), '\n')
  
  # Different true effects
  N.seq <- c(200, 400, 600, 800)
  models <- paste0("N=", N.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")
  
  
  res.alt <- NULL
  cc_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                     dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B <- matrix(rep(NA, K*length(models)), ncol = length(models),
                dimnames = list(NULL, models))
    #beta.n <- as.matrix(rnorm(K, sd = sqrt(var.b)))
    beta.n <- rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
    B[, models] <- beta.n %*% sqrt(N.seq)
    
    for(model in models){
      Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
      
      # Oracle
      method.tmp <- "Oracle"
      T.oracle <- as.numeric(Z.alt %*% sigma_trunc_inv %*% B[, model])
      cc_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% sigma_trunc_inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen_vec[, 1]
      cc_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_vec[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      cc_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      cc_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(cc_alt_all, file_p_alt)}
  }
  
}else if(change == "caus"){
  cat(paste0("Various ", change), '\n')
  
  caus.seq <- c(1, 5, 10, 30, 50)/100
  caus.num.seq <- floor(caus.seq * K)
  models <- paste0("caus=", caus.seq)
  
  res.alt <- NULL
  cc_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                     dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B <- matrix(rep(NA, K*length(models)), ncol = length(models),
                dimnames = list(NULL, models))
    B[, models] <- sapply(caus.num.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b)), rep(0, K-x))) ) * sqrt(N)
    
    for(model in models){
      Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
      
      # Oracle
      method.tmp <- "Oracle"
      T.oracle <- as.numeric(Z.alt %*% sigma_trunc_inv %*% B[, model])
      cc_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% sigma_trunc_inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen_vec[, 1]
      cc_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_lamb[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      cc_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      cc_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(cc_alt_all, file_p_alt)}
  }
  
}else if(change == "var"){
  cat(paste0("Various ", change), '\n')
  
  # Different true effects
  var.seq <- c(0.002, 0.003, 0.004, 0.005, 0.006)
  models <- paste0("var=", var.seq)
  
  cc_alt_all <- array(dim = c(N.sample, 1, length(models), N.sim),
                     dimnames = list(NULL, c("archie"), models, NULL))
  selected_gene_alt_all <- array(dim = c(N.sample, K, length(models), N.sim),
                                 dimnames = list(NULL, NULL, models, NULL))
  for(i in 1:N.sim){
    beta.n <- matrix(rep(NA, K*length(models)), ncol = length(models),
                dimnames = list(NULL, models))
    beta.n[, models] <- sapply(var.seq, function(x)
      as.matrix(c(rnorm(floor(caus*K), sd = sqrt(x)), rep(0, K-floor(caus*K)) ) ) )
    
    for(model in models){
      # archie
      Sigma_EE <- as.matrix(1)
      Sigma_GG <- Sigma
      Sigma_GE <- beta.n[, model, drop = FALSE] / sqrt(as.numeric(str_extract(model, '0\\.\\d+')))
      
      
      # us, vs, qs, K
      # us and vs in the output represent the selected SNPs mapped to corresponding selected genes (1 indiciating selected, 0 indicating not selected). 
      # In the output, note that SNPs 1 to 5 has been correctly mapped to genes 1 to 10 and SNPs 11 to 15 mapped to genes 11 to 20. 
      # ARCHIE usually numerically determines the number of sparse canonical correlation components to be extracted, however that is customizable as well. 
      # q is the aggregated measure of association explained by the corresponding component.
      res.alt <- archie_work(Sigma_GE, Sigma_GG, Sigma_EE, verbose = FALSE)
      
      cc_alt_all[, 'archie', model, i] <- max(res.alt$qs)
      selected_gene_alt_all[, , model, i] <- res.alt$us[, which.max(res.alt$qs)]
      
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(cc_alt_all, file_cc_alt)
      saveRDS(selected_gene_alt_all, file_selected_gene_alt)
    }
  }
  
}else{cat("Wrong parameter.")}


# print out key message or write out -----
saveRDS(cc_alt_all, file_cc_alt)
saveRDS(selected_gene_alt_all, file_selected_gene_alt)

