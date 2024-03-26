##############################################
########### fisher in simulation ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)


# I/O & paras -----
file_Sigma <- '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds'
K_all_gene <- 12102

p_cutoff_seq <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.05)

caus.seq <- c(1, 5, 10, 30, 50)/100
var.b <- 0.001
N <- 500
N.sample <- 1
N.sim <- 1e+5


## output -----
file_p_enrich <- str_glue('p_enrich_varb{var.b}_N{N}.rds')
file_p_plot <- str_glue('p_plot_varb{var.b}_N{N}.pdf')


# read files -----
Sigma <- as.matrix(readRDS(file_Sigma))
K <- dim(Sigma)[1]
K_null <- K_all_gene - K


# various caus as models -----
caus.num.seq <- floor(caus.seq * K)
models <- paste0("caus=", caus.seq)


# store enrich p for N.sample snps, various caus as models, and N.sim simulations
p_enrich <- array(
  dim = c(N.sample, length(p_cutoff_seq), length(models), N.sim),
  dimnames = list(NULL, str_glue('p_{p_cutoff_seq}'), models, NULL)
)

for(i in 1:N.sim){
  # simulate effect size b across various caus
  B <- matrix(rep(NA, K*length(models)), ncol = length(models),
              dimnames = list(NULL, models))
  B[, models] <- sapply(caus.num.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b)), rep(0, K-x))) ) * sqrt(N)
  
  
  for(model in models){
    # Z.alt for z of gene module with caus gene, using Sigma as module correlation
    Z.alt <- mvtnorm::rmvnorm(N.sample, B[, model], Sigma)
    
    
    # Z.null for z of all the other genes with non caus gene, genes are independent
    # as it's impossbile to simulation correlation for these non-caus genes
    Z.null <- rnorm(K_null*N.sample, 0, 1) %>% matrix(ncol = K_null)
    
    
    # convert z to p
    p.alt <- pchisq(Z.alt^2, df = 1, lower.tail = FALSE)
    p.null <- pchisq(Z.null^2, df = 1, lower.tail = FALSE)
    
    
    
    # fisher input
    aa <- lapply(p_cutoff_seq, function(x) {
      rowSums(p.alt < x)
    }) %>%
      setNames(str_glue('p_{p_cutoff_seq}')) %>%
      bind_cols()
    
    bb <- K - aa
    
    cc <- lapply(p_cutoff_seq, function(x) {
      rowSums(p.null < x)
    }) %>%
      setNames(str_glue('p_{p_cutoff_seq}')) %>%
      bind_cols()
    
    dd <- K_null - bb
    
    
    # convert fisher input to hypergeometric test input
    # as hypergeometric allows simulataneous calculation of N.sample tests
    q = aa
    m = aa + cc
    n = bb + dd
    k = aa + bb
    
    
    # store hypergeometric test for all N.sample and i-th simulation
    p_enrich[, str_glue('p_{p_cutoff_seq}'), model, i] <- lapply(str_glue('p_{p_cutoff_seq}'), function(x) {
      phyper(q[[x]], m[[x]], n[[x]], k[[x]], lower.tail = FALSE, log.p = FALSE)
    }) %>%
      setNames(str_glue('p_{p_cutoff_seq}')) %>%
      bind_cols() %>%
      unlist()
    
  }
  if(i %% 100 == 0) cat("Simulation: ", i, "\n")
}

saveRDS(p_enrich, file_p_enrich)




# plot of enrich p -----
plt_dat <- as.data.frame.table(p_enrich) %>%
  mutate(Var1 = NULL, Var4 = NULL) %>%
  rename('p_cutoff' = 'Var2', "model" = 'Var3', 'p_enrich' = 'Freq')

ggplot(plt_dat) +
  geom_boxplot(aes(x = model, y = -log10(p_enrich), color = p_cutoff))

# save
ggsave(file_p_plot, width = 6, height = 3)



# # qq plot & histogram & power -----
# n_test <- N.sim
# 
# # order statistic of null p
# expected <- seq(1, n_test) / (n_test+1)
# lexp <- -log10(expected)
# 
# ci_l <- -log10( qbeta(p = (1 - ci_level) / 2, shape1 = 1:n_test, shape2 = n_test:1) )
# ci_r <- -log10( qbeta(p = (1 + ci_level) / 2, shape1 = 1:n_test, shape2 = n_test:1) )
# 
# # obs
# observed <- apply(p_enrich, c(2, 3), sort) %>% as.data.frame.table()
# lobs <- -log10(observed)
# 
# 
# # data for plt
# df_plt = cbind(data.frame(x = lexp, ci_l = ci_l, ci_r = ci_r), lobs) %>%
#   pivot_longer(-c(x, ci_l, ci_r), names_to = "Type", values_to = "y")
# 
# 
# # qq plot
# ggplot(df_plt, aes(x = x, y = y, group = Type)) +
#   geom_ribbon(aes(ymin = ci_l, ymax = ci_r), fill = "#e5e5e5", color = "#e5e5e5") +
#   geom_abline(slope = 1, intercept = 0, color = "#595959", size = 0.7) +
#   geom_point(aes(color = Type), size = 0.5) +
#   labs(
#     x = bquote(Expected -log[10]~italic((P))),
#     y = bquote(Observed -log[10]~italic((P)))
#   ) +
#   scale_color_manual(
#     values = c(
#       "#85192d", "#0028a1", "#e89c31",
#       RColorBrewer::brewer.pal(8, "Dark2"),
#       RColorBrewer::brewer.pal(8, "Set1")
#     ) 
#   )

