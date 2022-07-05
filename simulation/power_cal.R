###############################################################
########### Calculate power using alt p and null p by FDR ###########
###############################################################
rm(list = ls())
library(data.table)

parameter <- commandArgs(trailingOnly = T)

file_p_alt <- parameter[1] #'new_Sigma/simulation.alt.N.lambda0.1.varb1e-3.K101.rds' # "simulation.alt.N.lambda0.1.K100.rds"
file_p_null <- parameter[2] #'new_Sigma/simulation.null.lambda0.1.K101.rds' # "simulation.null.lambda0.1.K100.rds"
file_out <- parameter[3] #'new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds' # "power.N.lambda0.1.K100.rds"
fdr_level <- as.numeric(parameter[4])
is_plot <- parameter[5]


p_alt_all <- readRDS(file_p_alt)
p_null_all <- readRDS(file_p_null)


N.sample <- dim(p_alt_all)[1]
C <- length(p_null_all[[1]])/N.sample
res <- NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  p.null <- data.table(paste0("null", 1:length(p_null_all[[paste0("p.null.", method.tmp)]])), p_null_all[[paste0("p.null.", method.tmp)]])
  res[[method.tmp]] <- apply(p_alt_all[ , method.tmp, , , drop = FALSE], c(3, 4), function(x) {
    p.obs = data.table(paste0("obs", 1:N.sample), as.numeric(x) )
    p.obs.rank = frank(p.obs, V2)
    names(p.obs.rank) = p.obs$V1

    all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
    names(all.rank) = c(p.obs$V1, p.null$V1)
    q = pmin((all.rank[names(p.obs.rank)]/p.obs.rank-1)/C, 1)
    
    sum(q < fdr_level)/N.sample
  } )
  cat("Method ", method.tmp, "is done.", '\n')
  saveRDS(res, file_out)
}

