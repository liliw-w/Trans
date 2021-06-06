rm(list = ls())
parameter = commandArgs(trailingOnly = T)
file_p_alt = parameter[1] # "simulation.alt.N.lambda0.1.K100.rds"
file_p_null = parameter[2] # "simulation.null.lambda0.1.K100.rds"
file_out = parameter[3] # "power.N.lambda0.1.K100.rds"

require(data.table)
p_alt_all = readRDS(file_p_alt)
p_null_all = readRDS(file_p_null)

N.sample = dim(p_alt_all)[1]
C = length(p_null_all[[1]])/N.sample
res = NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  p.null = data.table(paste0("null", 1:length(p_null_all[[paste0("p.null.", method.tmp)]])), p_null_all[[paste0("p.null.", method.tmp)]])
  res[[method.tmp]] = apply(p_alt_all[ , method.tmp, , 1:900], c(2, 3), function(x) {
    p.obs = data.table(paste0("obs", 1:N.sample), as.numeric(x) )
    p.obs.rank = frank(p.obs, V2)
    names(p.obs.rank) = p.obs$V1

    all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
    names(all.rank) = c(p.obs$V1, p.null$V1)
    q = pmin((all.rank[names(p.obs.rank)]/p.obs.rank-1)/C, 1)

    sum(q < 0.05)/N.sample
  } )
  cat("Method ", method.tmp, "is done.", '\n')
  saveRDS(res, file_out)

}
