###############################################################
########### Calculate power using alt p and null p by FDR ###########
###############################################################
# load packages -----
library(data.table)
library(tidyverse)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    /scratch/midway3/liliw1/paper1_sim/simulation_alt_lowCaus_highb_changecaus_varb0.2_N500_K101.rds
    /project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/simulation.null.lambda0.1.K101.rds
    0.05
    /scratch/midway3/liliw1/paper1_sim/power_lowCaus_highb_changecaus_varb0.2_N500_K101.rds
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_p_alt <- args[1] #'new_Sigma/simulation.alt.N.lambda0.1.varb1e-3.K101.rds' # "simulation.alt.N.lambda0.1.K100.rds"
file_p_null <- args[2] #'new_Sigma/simulation.null.lambda0.1.K101.rds' # "simulation.null.lambda0.1.K100.rds"
fdr_level <- as.numeric(args[3])

## output -----
file_out <- args[4] #'new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds' # "power.N.lambda0.1.K100.rds"


# read files -----
p_alt_all <- readRDS(file_p_alt)
p_null_all <- readRDS(file_p_null)


# power using empirical FDR correction -----
N.sample <- dim(p_alt_all)[1]
C <- length(p_null_all[[1]])/N.sample
res <- NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  p.null <- data.table(paste0("null", 1:length(p_null_all[[paste0("p.null.", method.tmp)]])), p_null_all[[paste0("p.null.", method.tmp)]])
  
  MyArray <- p_alt_all[ , method.tmp, , , drop = FALSE]
  res[[method.tmp]] <- apply(MyArray, 3, function(x) {
    tmp_p = apply(x, 3, function(y) list(as.numeric(y)))
    
    pbmcapply::pbmclapply(tmp_p, function(z) {
      p.obs = data.table(paste0("obs", 1:N.sample), unlist(z) )
      p.obs.rank = frank(p.obs, V2)
      names(p.obs.rank) = p.obs$V1
      
      all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
      names(all.rank) = c(p.obs$V1, p.null$V1)
      q = pmin((all.rank[names(p.obs.rank)]/p.obs.rank-1)/C, 1)
      
      sum(q < fdr_level)/N.sample
    }
    , mc.cores = 20) %>%
      list_c() %>%
      list()
  }) %>%
    setNames(dimnames(MyArray)[[3]])
  
  
  # print out key message or write out -----
  cat("Method ", method.tmp, "is done.", '\n')
  saveRDS(res, file_out)
}

