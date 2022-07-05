###############################################################
########### Simulation where beta and u1 align in direction ###########
########### PC1 has the best power ###########
###############################################################
rm(list = ls())
library(mvtnorm)

change <- 'N'
oracle.thre <- 0.1
PCO.script <- '/home/liliw1/Trans/simulation/script_lambda0.1/'
file.Sigma <- '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds'
file.res <- 'simulation.alt.pc1.rds'


var.b <- 0.001
caus <- 0.3
N <- 500

source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))


N.sample <- 10^4; N.sim <- 100


Sigma <- as.matrix(readRDS(file.Sigma))
SigmaO <- ModifiedSigmaOEstimate(Sigma)
K <- dim(Sigma)[1]
eigen.res <- eigen(Sigma)
lambdas <- eigen.res$values
eigen.vec <- eigen.res$vectors
Sigma.trunc.inv <- eigen.vec[, which(lambdas>oracle.thre)] %*% diag(1/lambdas[lambdas>oracle.thre]) %*% t(eigen.vec[, which(lambdas>oracle.thre)])

if(change == "N"){
  cat(paste0("Various ", change), '\n')
  
  # Different true effects
  N.seq <- c(200, 400, 600, 800)
  models <- paste0("N=", N.seq)
  
  
  res.alt <- NULL
  p_alt_all <- array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B <- matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    beta.n <- as.matrix(rnorm(K, sd = sqrt(var.b)))
    beta.n <- rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
    B[, models] <- beta.n %*% sqrt(N.seq)
    B[, models] <- eigen.vec[, 1] %*% t( as.matrix(sqrt(colSums(B[, models]^2)) ))/caus*1.2
    
    for(model in models){
      Z.alt <- rmvnorm(N.sample, B[, model], Sigma)
      
      # Oracle
      method.tmp <- "Oracle"
      T.oracle <- as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] <- 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)
      
      # PC1
      method.tmp <- "PC1"
      PC1 <- Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] <- as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
      
      # univariate minp
      method.tmp <- "minp"
      p_alt_all[, method.tmp, model, i] <- apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      
      # PCO
      method.tmp <- "PCO"
      p_alt_all[, method.tmp, model, i] <- as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file.res)}
  }
  
}

saveRDS(p_alt_all, file.res)
