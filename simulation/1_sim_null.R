##############################################
########### Run null simulations to see association methods are well calibrated ###########
##############################################
# load packages -----
rm(list = ls())
library(mvtnorm)


# I/O & paras -----
n_sim <- 10^7
dir_pco <- './script_lambda0.1/'
file_Sigma <- '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds'

## output -----
file_dat_null <- "/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/simulation.null.lambda0.1.K101.rds"


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



# simulate null z-scores -----
p_null_all  <- list()
z_null <- rmvnorm(n_sim, rep(0, K), Sigma)


# prepare eigenvalues and eigenvectors as input for three methods -----
SigmaO <- ModifiedSigmaOEstimate(Sigma)
eigen_res <- eigen(Sigma)
eigen_lamb <- eigen_res$values
eigen_vec <- eigen_res$vectors


# run three association test and calculate p-value -----
# PCO
p_null_all$'p.null.PCO' <- ModifiedPCOMerged(
  Z.mat = z_null, Sigma = Sigma, SigmaO = SigmaO
) |> as.numeric()

cat("PCO done. \n\n")


# PC1
PC1 <- z_null %*% eigen_vec[, 1]
p_null_all$'p.null.PC1' <- 2*pnorm(-abs(PC1/sqrt(eigen_lamb[1])))|> as.numeric()

cat("PC1 done. \n\n")


# univariate minp
p_null_all$'p.null.minp' <- apply(z_null, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

cat("minp done. \n\n")


# print out key message or write out -----
saveRDS(p_null_all, file = file_dat_null)

