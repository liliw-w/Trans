###############################################################
########### Power plot when low caus, high var ###########
########### for all varb ###########
###############################################################
# load packages -----
rm(list = ls())
library(tidyverse)
library(ggpubr)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '

    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_power_all <- list.files(
  '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma',
  'power_lowCaus_highb_changecaus_varb0[.]\\d+_N500_K101.rds',
  full.names = TRUE
)
file_fig_power <- '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/plt_power_lowCaus_highb_changecaus_allvarb_N500_K101.pdf'


## output -----


# read & organize files -----
res.alt <- lapply(
  file_power_all,
  function(x) {
    map_dfr(
      readRDS(x),
      ~bind_cols(.x) %>% 
        setNames(names(.x)) %>%
        pivot_longer(everything(), names_to = "model", values_to = "power") %>%
        mutate("model" = as.numeric(str_extract(model, "\\d+.*\\d+$"))),
      .id = "method"
    )
  }
) %>% 
  setNames(str_extract(file_power_all, "varb0[.]\\d+") %>% str_extract('0[.]\\d+')) %>% 
  bind_rows(.id = "varb")

res.alt$method <- factor(res.alt$method, c("PCO", "PC1", "minp"), c("Trans-PCO", "PC1", "MinP"))
res.alt$model <- factor(res.alt$model, levels = sort(unique(res.alt$model)))
res.alt$varb <- factor(
  res.alt$varb, 
  levels = sort(unique(as.numeric(res.alt$varb))),
  labels = paste0('varb=', sort(unique(as.numeric(res.alt$varb))))
)


# print out key message or write out -----
plt_dat <- filter(res.alt, varb != "varb=0.005") %>%
  filter(method != "PC1")

ggline(
  plt_dat, 
  x = "model", y = "power", add = "mean_se",
  facet.by = "varb",
  color = "method", palette = "jco"
) +
  stat_compare_means(
    aes(group = method), label = "p.format", 
    label.y = c(0.1, 0.5, 0.7)
  ) +
  labs(
    x = "Causal Proportion", 
    y = "Power", 
    color = "Method"
  )


# compare_means(power ~ method, data = plt_dat, group.by = c("model", "varb")) %>%
#   arrange(varb, model)

# save figure -----
ggsave(
  file_fig_power, 
  height = 4, width = 5
)

