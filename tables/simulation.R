##############################################
########### write out simulation numerical results ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)


# I/O & paras -----
file_power_all_seq <- c(
  'new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds',
  'new_Sigma/power.caus.lambda0.1.varb1e-3.K101.rds',
  'new_Sigma/power.var.lambda0.1.varb1e-3.K101.rds'
)
model_name_seq <- c(
  "Sample Size",
  "Causal Proportion",
  "Genetic Variance"
)

## output -----
file_simulation_res <- 'new_Sigma/numerical_res.txt'


# iterate through simulations and get numerical results -----
res <- NULL
for(s in 1:length(file_power_all_seq)){
  
  # read files -----
  file_power_all <- file_power_all_seq[s]
  model_name <- model_name_seq[s]
  
  res.alt <- map_dfr(
    readRDS(file_power_all),
    ~as_tibble(.x, rownames = "paras") %>%
      pivot_longer(!paras, names_to = NULL, values_to = "power") %>%
      mutate("paras" = str_extract(paras, "\\d+.*\\d+$")),
    .id = "method"
  )
  res.alt$method <- factor(res.alt$method, c("PCO", "PC1", "minp"), c("Trans-PCO", "PC1", "MinP"))
  
  
  # calcaulate mean and ci -----
  res <- rbind(
    res,
    group_by(res.alt, paras, method) %>%
      summarise(
        mean_cl_normal(power, conf.int = .95) %>%
          rename(power_mean = y, lower_ci = ymin, upper_ci = ymax)
      ) %>%
      ungroup() %>%
      mutate("scenario" = model_name) %>%
      relocate(scenario, .before = everything())
  )
}

# print out key message or write out -----
fwrite(
  res,
  file_simulation_res,
  quote = FALSE, sep = "\t"
)

