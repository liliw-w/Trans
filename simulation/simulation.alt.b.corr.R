# Z alternative
# FDR based on empirical pnull

rm(list = ls())
parameter = commandArgs(trailingOnly = T)

change = 'N'
oracle.thre = 0.1
PCO.script = "/home/liliw1/Trans/simulation/script_lambda0.1/"
file.Sigma = "~/xuanyao_llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds" #'Sigma.DGN.module13_chr3.100.rds'
file.res = 'simulation.alt.N.lambda0.1.K100_b_indep_caus1.rds'

var.b = 0.001
caus = 1 # 0.3
N = 500

library(mvtnorm)

source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))


N.sample = 10^4; N.sim = 10^3
#thre = 0.05/(50*1000000)

Sigma = as.matrix(readRDS(file.Sigma))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors
Sigma.trunc.inv = eigen.vec[, which(lambdas>oracle.thre)] %*% diag(1/lambdas[lambdas>oracle.thre]) %*% t(eigen.vec[, which(lambdas>oracle.thre)])

#if(change == "N"){
  cat(paste0("Various ", change), '\n')
  
  # Different true effects
  N.seq = c(200, 400, 600, 800)
  models = paste0("N=", N.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")
  
  
  res.alt = NULL
  p_alt_all = array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    #beta.n = as.matrix(rnorm(K, sd = sqrt(var.b)))
    
    #var.b.seq <- matrix(c(rep(var.b, floor(caus*K)), rep(0, K-floor(caus*K)) ))
    #var.b.mat <- sqrt(var.b.seq) %*% t(sqrt(var.b.seq))
    #beta.n <- rmvnorm(1, rep(0, K), Sigma * var.b.mat)
    #beta.n <- t(beta.n)
    
    beta.n = rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))),
                   as.matrix(rep(0, K-floor(caus*K))))
    B[, models] = beta.n %*% sqrt(N.seq)
    
    
    for(model in models){
      #model <- models[4]
      
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)
      
      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] = 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)
      #sum(p_alt_all[, method.tmp, model, i] < 1e-6)/N.sample
      
      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
      #sum(p_alt_all[, method.tmp, model, i] < 1e-5)/N.sample
      
      # univariate minp
      method.tmp = "minp"
      p_alt_all[, method.tmp, model, i] = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      #sum(p_alt_all[, method.tmp, model, i] < 1e-5)/N.sample
      
      # PCO
      #method.tmp = "PCO"
      #p_alt_all[, method.tmp, model, i] = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      #sum(p_alt_all[, method.tmp, model, i] < 1e-6)/N.sample
      
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file.res)}
  }
  
  
  #apply(p_alt_all[, , , 1], c(2, 3), function(x) sum(x<0.05))
  #apply(p_alt_all_old[, , , 1], c(2, 3), function(x) sum(x<0.05))
  
  #p_alt_all_old <- readRDS('~/xuanyao_llw/simulation_lambda0.1/new_Sigma/simulation.alt.N.lambda0.1.varb1e-3.K101.rds')
  
  
  
#}
saveRDS(p_alt_all, file.res)
