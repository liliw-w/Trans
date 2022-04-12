### Supp table of T1E in simulations
rm(list = ls())
file_null_p = 'simulation.null.lambda0.1.K101.rds'
alpha_seq = c(0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6)
file_T1E = "T1E_table.txt"

null_p = readRDS(file_null_p)
str(null_p)

f_compare <- function(x, y, null_p, alpha_seq){
  sum(null_p[[x]] < y)/length(x)
}

T1E = sapply(alpha_seq, function(x) colMeans(as.matrix(as.data.frame(null_p)) < x) )
colnames(T1E) = alpha_seq

write.table(format(T1E, digits = 2),
            file = file_T1E,
            quote = FALSE,
            sep = "\t")

