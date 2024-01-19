##############################################
########### fisher power in simulation ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')

theme_set(
  theme_my_pub(legend.position = "bottom")
)


# I/O & paras -----
file_Sigma <- '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds'
K_all_gene <- 12102
fdr_level <- 0.1

p_cutoff_seq <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.05)

caus.seq <- c(1, 5, 10, 30, 50)/100
var.b <- 0.001
N <- 500
N.sample <- 1e+4
N.sim <- 1e+3



## output -----
file_power_enrich <- str_glue('power_enrich_varb{var.b}_N{N}.rds')
file_aa_all <- str_glue('aa_varb{var.b}_N{N}.rds')
file_cc_all <- str_glue('cc_varb{var.b}_N{N}.rds')
file_power_plot <- str_glue('power_plot_varb{var.b}_N{N}.pdf')
file_contigency_table <- str_glue('table_enrich_varb{var.b}_N{N}.pdf')


# read files -----
Sigma <- as.matrix(readRDS(file_Sigma))
K <- dim(Sigma)[1]
K_null <- K_all_gene - K


# various caus as models -----
caus.num.seq <- floor(caus.seq * K)
models <- paste0("caus=", caus.seq)


# store enrich power and contigency table for N.sample snps, various caus as models, and N.sim simulations
power_enrich <- array(
  dim = c(length(p_cutoff_seq), length(models), N.sim),
  dimnames = list(str_glue('p_{p_cutoff_seq}'), models, NULL)
)
aa_all <- array(
  dim = c(length(p_cutoff_seq), length(models), N.sim),
  dimnames = list(str_glue('p_{p_cutoff_seq}'), models, NULL)
)
cc_all <- array(
  dim = c(length(p_cutoff_seq), length(models), N.sim),
  dimnames = list(str_glue('p_{p_cutoff_seq}'), models, NULL)
)

for(i in 1:N.sim){
  # simulate effect size b across various caus
  B <- matrix(rep(NA, K*length(models)), ncol = length(models),
              dimnames = list(NULL, models))
  B[, models] <- sapply(caus.num.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b)), rep(0, K-x))) ) * sqrt(N)
  
  
  # Z.null for z of all the other genes with non caus gene, genes are independent
  # as it's impossbile to simulation correlation for these non-caus genes
  Z.null <- rnorm(K_null*N.sample, 0, 1) %>% matrix(ncol = K_null)
  
  # convert z to p
  p.null <- pchisq(Z.null^2, df = 1, lower.tail = FALSE)
  
  
  for(model in models){
    # Z.alt for z of gene module with caus gene, using Sigma as module correlation
    Z.alt <- mvtnorm::rmvnorm(N.sample, B[, model], Sigma)
    
    
    # convert z to p
    p.alt <- pchisq(Z.alt^2, df = 1, lower.tail = FALSE)
    
    
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
    
    dd <- K_null - cc
    
    
    # convert fisher input to hypergeometric test input
    # as hypergeometric allows simulataneous calculation of N.sample tests
    q = aa
    m = aa + cc
    n = bb + dd
    k = aa + bb
    
    
    # store hypergeometric test for all N.sample and i-th simulation
    p_enrich <- lapply(str_glue('p_{p_cutoff_seq}'), function(x) {
      phyper(q[[x]], m[[x]], n[[x]], k[[x]], lower.tail = FALSE, log.p = FALSE)
    }) %>%
      setNames(str_glue('p_{p_cutoff_seq}')) %>%
      bind_cols()
    
    # manually assign those with aa == 0
    ind_0 <- aa == 0
    p_enrich[ind_0] <- 1
    
    
    # qvalue & power
    power_enrich[str_glue('p_{p_cutoff_seq}'), model, i] <- apply(
      p_enrich, 
      2, 
      function(x) {
        sum(qvalue::qvalue(x, fdr.level = fdr_level, lambda = 0)$significant)/N.sample
      }
    )
    
    # fisher contigency table
    aa_all[str_glue('p_{p_cutoff_seq}'), model, i] <- colMeans(aa)
    cc_all[str_glue('p_{p_cutoff_seq}'), model, i] <- colMeans(cc)
    
  }
  
  cat("Simulation: ", i, "\n")
  
  if(i %% 20 == 0) {
    saveRDS(power_enrich, file_power_enrich)
    saveRDS(aa_all, file_aa_all)
    saveRDS(cc_all, file_cc_all)
  }
}


# plot of enrich power -----
plt_dat <- as.data.frame.table(power_enrich) %>%
  mutate(Var3 = NULL) %>%
  rename('p_cutoff' = 'Var1', "model" = 'Var2', 'power_enrich' = 'Freq')

ggpubr::ggarrange(
  ggplot(plt_dat) +
    geom_boxplot(
      aes(x = model, y = power_enrich, color = p_cutoff),
      outlier.alpha = 0.2, outlier.size = 0.1
    ) +
    scale_color_brewer(palette = "Dark2"),
  
  ggplot(plt_dat) +
    facet_wrap(~model, nrow = 2) +
    geom_freqpoly(
      aes(x = power_enrich, color = p_cutoff),
    ) +
    scale_color_brewer(palette = "Dark2"),
  
  nrow = 2, heights = c(1, 1.5)
)


# save
ggsave(file_power_plot, width = 6, height = 6)



# plot of fisher contigency table -----
# aa_all <- readRDS(file_aa_all)
# cc_all <- readRDS(file_cc_all)


plt_dat_aa <- as.data.frame.table(aa_all) %>%
  mutate(Var3 = NULL) %>%
  rename('p_cutoff' = 'Var1', "model" = 'Var2', 'power_enrich' = 'Freq')
plt_dat_cc <- as.data.frame.table(cc_all) %>%
  mutate(Var3 = NULL) %>%
  rename('p_cutoff' = 'Var1', "model" = 'Var2', 'power_enrich' = 'Freq')

ggpubr::ggarrange(
  ggplot(plt_dat_aa) +
    geom_boxplot(
      aes(x = model, y = power_enrich, color = p_cutoff),
      outlier.alpha = 0.2, outlier.size = 0.1
    ) +
    labs(y = "signals in module (a)") +
    scale_color_brewer(palette = "Dark2"),
  
  ggplot(plt_dat_cc) +
    geom_boxplot(
      aes(x = model, y = power_enrich, color = p_cutoff),
      outlier.alpha = 0.2, outlier.size = 0.1
    ) +
    labs(y = "signals outside module (c)") +
    scale_color_brewer(palette = "Dark2"),
  
  nrow = 2, heights = c(1, 1)
)

ggsave(file_contigency_table, width = 6, height = 4)


