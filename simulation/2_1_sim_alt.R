##############################################
############ p-values under alternative ############
##############################################
# load packages -----
rm(list = ls())
library(mvtnorm)


# I/O & paras -----
parameter <- commandArgs(trailingOnly = T)

change <- parameter[1] #'N'
oracle_thre <- as.numeric(parameter[2]) #0.1
dir_pco <- parameter[3] #'./script_lambda0.1/'
file_Sigma <- parameter[4] #'Sigma.DGN.module13_chr3.100.rds'
file_p_alt <- parameter[5] #'simulation.alt.N.lambda0.1.K100.rds'

var.b <- as.numeric(parameter[6]) #0.001
caus <- as.numeric(parameter[7]) #0.3
N <- as.numeric(parameter[8]) #500

N.sample <- as.numeric(parameter[9]) #10^4
N.sim <- as.numeric(parameter[10]) #10^3



# read files -----
# source pco test
source(paste0(dir_pco, "ModifiedPCOMerged.R"))
source(paste0(dir_pco, "liu.R"))
source(paste0(dir_pco, "liumod.R"))
source(paste0(dir_pco, "davies.R"))
dyn.load(paste0(dir_pco, "qfc.so"))
source(paste0(dir_pco, "ModifiedSigmaOEstimate.R"))

Sigma <- as.matrix(readRDS(file_Sigma))
K <- dim(Sigma)[1]


# prepare eigenvalues and eigenvectors as input for three methods -----
SigmaO <- ModifiedSigmaOEstimate(Sigma)
eigen_res <- eigen(Sigma)
eigen_lamb <- eigen_res$values
eigen_vec <- eigen_res$vectors


sigma_trunc_inv <- eigen_vec[, which(eigen_lamb>oracle_thre)] %*% diag(1/eigen_lamb[eigen_lamb>oracle_thre]) %*% t(eigen_vec[, which(eigen_lamb>oracle_thre)])


# p values of alt z's under various N -----
if(change == "N"){
  cat(paste0("Various ", change), '\n')
  
  # Different true effects
  N.seq <- c(200, 400, 600, 800)
  models <- paste0("N=", N.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")
  
  
  res.alt <- NULL
  p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
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
      p_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% sigma_trunc_inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen_vec[, 1]
      p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_vec[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      p_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file_p_alt)}
  }
  
}else if(change == "caus"){
  cat(paste0("Various ", change), '\n')
  
  caus.seq <- c(1, 5, 10, 30, 50)/100
  caus.num.seq <- floor(caus.seq * K)
  models <- paste0("caus=", caus.seq)
  
  res.alt <- NULL
  p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
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
      p_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% sigma_trunc_inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen_vec[, 1]
      p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_lamb[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      p_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file_p_alt)}
  }
  
}else if(change == "caus_fix"){
  cat(paste0("Various ", change), '\n')
  
  caus.seq <- c(1, 5, 10, 30, 50)/100
  caus.num.seq <- floor(caus.seq * K)
  models <- paste0("caus.fix=", caus.seq)
  
  res.alt <- NULL
  p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                     dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B <- matrix(rep(NA, K*length(models)), ncol = length(models),
                dimnames = list(NULL, models))
    B[, models] <- sapply(caus.num.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b/x)), rep(0, K-x))) ) * sqrt(N)
    
    for(model in models){
      Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
      
      # Oracle
      method.tmp <- "Oracle"
      T.oracle <- as.numeric(Z.alt %*% sigma_trunc_inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% sigma_trunc_inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen_vec[, 1]
      p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_lamb[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      p_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file_p_alt)}
  }
  
}else if(change == "var"){
  cat(paste0("Various ", change), '\n')
  
  # Different true effects
  var.seq <- c(0.002, 0.003, 0.004, 0.005, 0.006)
  models <- paste0("var=", var.seq)
  
  res.alt <- NULL
  p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                     dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B <- matrix(rep(NA, K*length(models)), ncol = length(models),
                dimnames = list(NULL, models))
    B[, models] <- sapply(var.seq, function(x)
      as.matrix(c(rnorm(floor(caus*K), sd = sqrt(x)), rep(0, K-floor(caus*K)) ) ) ) * sqrt(N)
    
    for(model in models){
      Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
      
      # Oracle
      method.tmp <- "Oracle"
      T.oracle <- as.numeric(Z.alt %*% sigma_trunc_inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% sigma_trunc_inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen_vec[, 1]
      p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(eigen_lamb[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      p_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file_p_alt)}
  }
  
}else{cat("Wrong parameter.")}


# print out key message or write out -----
saveRDS(p_alt_all, file_p_alt)
