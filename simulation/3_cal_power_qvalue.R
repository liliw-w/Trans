###############################################################
########### Calculate power using qvalue and fdr 10% ###########
###############################################################
# load packages -----
library(data.table)
library(tidyverse)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    /project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/simulation_alt_lowCaus_highb_changecaus_varb0.2_N500_K101.rds
    0.1
    /project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/power_qvalue_fdr10_lowCaus_highb_changecaus_varb0.2_N500_K101.rds
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_p_alt <- args[1] #'new_Sigma/simulation.alt.N.lambda0.1.varb1e-3.K101.rds' # "simulation.alt.N.lambda0.1.K100.rds"
fdr_level <- as.numeric(args[2])

## output -----
file_out <- args[3] #'new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds' # "power.N.lambda0.1.K100.rds"


# read files -----
p_alt_all <- readRDS(file_p_alt)


# power using empirical FDR correction -----
N.sample <- dim(p_alt_all)[1]
res <- NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  MyArray <- p_alt_all[ , method.tmp, , , drop = FALSE]
  
  
  res[[method.tmp]] <- apply(MyArray, 3, function(x) {
    tmp_p = apply(x, 3, function(y) list(as.numeric(y)))
    
    pbmcapply::pbmclapply(tmp_p, function(z) {
      sum(qvalue::qvalue(pmin(unlist(z), 1), fdr.level = fdr_level, lambda = 0)$significant)/N.sample
    }
    , mc.cores = 20
    ) %>%
      list_c() %>%
      list()
  }
  ) %>%
    setNames(dimnames(MyArray)[[3]])
  
  
  # print out key message or write out -----
  cat("Method ", method.tmp, "is done.", '\n')
  saveRDS(res, file_out)
}

