###############################################################
########### Simulation using low caus, high var ###########
########### to look at if uni outperforms multi ###########
###############################################################
# load packages -----
library(mvtnorm)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    caus
    0.1
    ./script_lambda0.1/
    ./new_Sigma/Sigma-new_DGN_module29_K101.rds
    ./new_Sigma/simulation_alt_lowCaus_highb_changecaus_varb0.2_N500_K101.rds
    0.2
    0.3
    500
    10000
    1000
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
file_p_alt <- args[5]

var.b <- as.numeric(args[6])
caus <- as.numeric(args[7])
N <- as.numeric(args[8])

N.sample <- as.numeric(args[9])
N.sim <- as.numeric(args[10])

# low caus (high sparsity)
caus.seq <- c(1, 5, 10)/100



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


# p values of alt z's under various models -----
cat(paste0("Various ", change), '\n')

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
  if(i %% 100 == 0){
    saveRDS(p_alt_all, file_p_alt)}
}

# print out key message or write out -----
saveRDS(p_alt_all, file_p_alt)
