##############################################
########### power plot of rotival et al. & other methods ###########
########### qvalue ###########
########### fdr 10% ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
library(ggpubr)


# I/O & paras -----
N <- 500

## power of other methods
file_power_other <- list.files(
  path = '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma',
  pattern = str_glue('power_qvalue_fdr10_.*.rds'),
  full.names = TRUE
)

## power of rotival et al
file_power_enrich <- list.files(
  path = ".",
  pattern = str_glue('power_enrich_varb\\d+[.]\\d+_N{N}.rds'),
  full.names = TRUE
)


## output -----
file_fig_power_rotival_others <- str_glue('fig_power_rotival_others_varball_N{N}.pdf')
file_fig_power_rotival_only <- str_glue('fig_power_rotival_only_varball_N{N}.pdf')



# read files -----
## power of other methods -----
power_other <- lapply(
  file_power_other,
  function(x) {
    map_dfr(
      readRDS(x),
      ~bind_cols(.x) %>% 
        setNames(names(.x)) %>%
        pivot_longer(everything(), names_to = "model", values_to = "power"),
      .id = "method"
    )
  }
) %>% 
  setNames(str_extract(file_power_other, "varb0[.]\\d+")) %>% 
  bind_rows(.id = "varb")



## power of rotival et al -----
power_enrich <- lapply(
  file_power_enrich,
  function(x) {
    readRDS(x) %>%
      as.data.frame.table() %>%
      mutate(Var3 = NULL) %>%
      rename('p_cutoff' = 'Var1', "model" = 'Var2', 'power' = 'Freq')
  }
) %>%
  setNames(str_extract(file_power_enrich, "varb\\d+[.]\\d+")) %>%
  bind_rows(.id = "varb") %>%
  mutate(method = str_glue("rotival_{p_cutoff}"))


## combined power -----
res <- bind_rows(power_other, select(power_enrich, !p_cutoff))




# data for plot -----
# factorize and reorder
res$model <- as.numeric(str_extract(res$model, "\\d+.*\\d+$"))
res$varb <- as.numeric(str_extract(res$varb, "\\d+.*\\d+$"))

res$model <- factor(res$model, levels = sort(unique(res$model)))
res$varb <- factor(
  res$varb, 
  levels = sort(unique(res$varb))
)


## power comparison plot of rotival et al & other methods -----
filter(res, method %in% c("PCO", "PC1", "minp", "rotival_p_1e-05")) %>%
  ggline(
    x = "model", y = "power", add = "mean_se",
    facet.by = "varb", ncol = 2,
    color = "method", palette = "jco"
  ) +
  labs(
    x = "Causal Proportion", 
    y = "Power", 
    color = "Method"
  ) +
  scale_colour_manual(
    breaks = c(
      "PCO", "PC1", "minp", 
      'rotival_p_1e-05', 'rotival_p_1e-04', 'rotival_p_0.001', 'rotival_p_0.01', 'rotival_p_0.05'
    ),
    values = c(
      "PCO" = "#85192d", "PC1" = "#1d349a", "minp" = "#e89c31", 
      'rotival_p_1e-05' = "#810F7C", 'rotival_p_1e-04' = "#8856A7", 'rotival_p_0.001' = "#8C96C6", 'rotival_p_0.01' = "#B3CDE3", 'rotival_p_0.05' = "#EDF8FB"
    )
  )


ggsave(file_fig_power_rotival_others, width = 5, height = 6)




## power comparison plot of rotival et al across all varb & p cutoff -----
filter(res, !(method %in% c("PCO", "PC1", "minp"))) %>%
  ggline(
    x = "model", y = "power", add = "mean_se",
    facet.by = "varb", ncol = 2,
    color = "method", palette = "jco"
  ) +
  labs(
    x = "Causal Proportion", 
    y = "Power", 
    color = "Method"
  ) +
  scale_colour_manual(
    breaks = c(
      "PCO", "PC1", "minp", 
      'rotival_p_1e-05', 'rotival_p_1e-04', 'rotival_p_0.001', 'rotival_p_0.01', 'rotival_p_0.05'
    ),
    values = c(
      "PCO" = "#85192d", "PC1" = "#1d349a", "minp" = "#e89c31", 
      'rotival_p_1e-05' = "#810F7C", 'rotival_p_1e-04' = "#8856A7", 'rotival_p_0.001' = "#8C96C6", 'rotival_p_0.01' = "#B3CDE3", 'rotival_p_0.05' = "#EDF8FB"
    )
  )

ggsave(file_fig_power_rotival_only, width = 5, height = 6)



# ## boxplot - power comparison plot of rotival et al across all varb & p cutoff -----
# ggpubr::ggarrange(
#   ggpubr::ggarrange(
#     filter(power_enrich, varb == 'varb0.001') %>%
#       ggplot() +
#       facet_wrap(~varb) +
#       geom_boxplot(
#         aes(x = model, y = power, color = p_cutoff),
#         outlier.alpha = 0.2, outlier.size = 0.1
#       ) +
#       labs(y = "Power") +
#       scale_color_brewer(palette = "Dark2"),
#     
#     ggplot() + theme_void(),
#     ncol = 2
#   ),
#   
#   
#   filter(power_enrich, varb != 'varb0.001') %>%
#     ggplot() +
#     facet_wrap(~varb) +
#     geom_boxplot(
#       aes(x = model, y = yy, color = p_cutoff),
#       outlier.alpha = 0.2, outlier.size = 0.1
#     ) +
#     labs(y = "Power") +
#     scale_color_brewer(palette = "Dark2"),
#   
#   ncol = 1, heights = c(1, 1.5), widths = c(1, 1.5)
# )
# 
# 
# file_power_plot <- str_glue('power_plot_varball_N{N}.pdf')
# ggsave(file_power_plot, width = 7, height = 7)
# 
# 
# 
# group_by(power_enrich, varb, p_cutoff, model) %>%
#   summarise('mean_power' = mean(yy, na.rm = TRUE))




# ## aa plot of rotival et al across all varb & p cutoff -----
# file_aa_all <- list.files(
#   pattern = str_glue('aa_varb\\d+[.]\\d+_N{N}.rds')
# )
# 
# aa_all <- lapply(
#   file_aa_all,
#   function(x) {
#     readRDS(x) %>%
#       as.data.frame.table() %>%
#       mutate(Var3 = NULL) %>%
#       rename('p_cutoff' = 'Var1', "model" = 'Var2', 'yy' = 'Freq')
#   }
# ) %>%
#   setNames(str_extract(file_aa_all, "varb\\d+[.]\\d+")) %>%
#   bind_rows(.id = "varb")
# 
# 
# ggpubr::ggarrange(
#   ggpubr::ggarrange(
#     filter(aa_all, varb == 'varb0.001') %>%
#       ggplot() +
#       facet_wrap(~varb) +
#       geom_boxplot(
#         aes(x = model, y = yy, color = p_cutoff),
#         outlier.alpha = 0.2, outlier.size = 0.1
#       ) +
#       labs(y = "signals in module (a)") +
#       scale_color_brewer(palette = "Dark2"),
#     
#     ggplot() + theme_void(),
#     ncol = 2
#   ),
#   
#   
#   filter(aa_all, varb != 'varb0.001') %>%
#     ggplot() +
#     facet_wrap(~varb) +
#     geom_boxplot(
#       aes(x = model, y = yy, color = p_cutoff),
#       outlier.alpha = 0.2, outlier.size = 0.1
#     ) +
#     labs(y = "signals in module (a)") +
#     scale_color_brewer(palette = "Dark2"),
#   
#   ncol = 1, heights = c(1, 1.5), widths = c(1, 1.5)
# )
# 
# 
# file_aa_plot <- str_glue('aa_varball_N{N}.pdf')
# ggsave(file_aa_plot, width = 7, height = 7)
# 
# 
# group_by(aa_all, varb, p_cutoff, model) %>%
#   summarise('mean_aa' = mean(yy, na.rm = TRUE))




# ## cc plot of rotival et al across all varb & p cutoff -----
# file_cc_all <- list.files(
#   pattern = str_glue('cc_varb\\d+[.]\\d+_N{N}.rds')
# )
# 
# cc_all <- lapply(
#   file_cc_all,
#   function(x) {
#     readRDS(x) %>%
#       as.data.frame.table() %>%
#       mutate(Var3 = NULL) %>%
#       rename('p_cutoff' = 'Var1', "model" = 'Var2', 'yy' = 'Freq')
#   }
# ) %>%
#   setNames(str_extract(file_cc_all, "varb\\d+[.]\\d+")) %>%
#   bind_rows(.id = "varb")
# 
# 
# 
# ggpubr::ggarrange(
#   ggpubr::ggarrange(
#     filter(cc_all, varb == 'varb0.001') %>%
#       ggplot() +
#       facet_wrap(~varb) +
#       geom_boxplot(
#         aes(x = model, y = yy, color = p_cutoff),
#         outlier.alpha = 0.2, outlier.size = 0.1
#       ) +
#       labs(y = "signals outside module (c)") +
#       scale_color_brewer(palette = "Dark2"),
#     
#     ggplot() + theme_void(),
#     ncol = 2
#   ),
#   
#   
#   filter(cc_all, varb != 'varb0.001') %>%
#     ggplot() +
#     facet_wrap(~varb) +
#     geom_boxplot(
#       aes(x = model, y = yy, color = p_cutoff),
#       outlier.alpha = 0.2, outlier.size = 0.1
#     ) +
#     labs(y = "signals outside module (c)") +
#     scale_color_brewer(palette = "Dark2"),
#   
#   ncol = 1, heights = c(1, 1.5), widths = c(1, 1.5)
# )
# 
# 
# file_cc_plot <- str_glue('cc_varball_N{N}.pdf')
# ggsave(file_cc_plot, width = 7, height = 7)
# 
# 
# group_by(cc_all, varb, p_cutoff, model) %>%
#   summarise('mean_cc' = mean(yy, na.rm = TRUE))

