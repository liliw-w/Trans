rm(list = ls())
library(tidyverse)


file_dat_null = "~/scratch/simulation_lambda0.1/new_Sigma/simulation.null.lambda0.1.K101.rds"

source("/home/liliw1/Trans/plot/theme_my_pub.R")

group_order = c("p.null.PCO", "p.null.PC1", "p.null.minp")
group_label = c("Trans-PCO", "PC1", "MinP")
ci_level = 0.95


# read data
p_null_all = readRDS(file_dat_null)



#reset obs p-values that are 0 to a fixed value, here I use the min non-zero p/10
input = as_tibble(p_null_all)
input[input == 0] = min(input[input != 0])/10

# number of samples
n = nrow(input)

# expected
expected = seq(1, n) / (n+1)
lexp = -log10(expected)

# order statistic of null p
ci_l = -log10( qbeta(p = (1 - ci_level) / 2, shape1 = 1:n, shape2 = n:1) )
ci_r = -log10( qbeta(p = (1 + ci_level) / 2, shape1 = 1:n, shape2 = n:1) )


# obs
observed = apply(input, 2, sort) %>% as.data.frame()
lobs = -log10(observed)


# take only a subset of null p's, to save image space
ind_sub = c(
  1:sum(lexp > 4),
  seq(from = sum(lexp > 4), to = sum(lexp > 2), length.out = 2000) %>% ceiling(),
  seq(from = sum(lexp > 2), to = n, length.out = 2000) %>% ceiling()
)


# data for plt
df_plt = cbind(data.frame(x = lexp, ci_l = ci_l, ci_r = ci_r), lobs) %>%
  slice(ind_sub) %>%
  pivot_longer(-c(x, ci_l, ci_r), names_to = "Type", values_to = "y")
# set group order in plt
group_order = if(is.null(group_order)) unique(df_plt$Type) else group_order
group_label = if(is.null(group_label)) group_order else group_label
df_plt$Type = factor(df_plt$Type, levels = group_order, labels = group_label)


# plt
base_plt <- ggplot(df_plt, aes(x = x, y = y, group = Type)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_r), fill = "#e5e5e5", color = "#e5e5e5") +
  geom_abline(slope = 1, intercept = 0, color = "#595959", size = 0.7) +
  geom_point(aes(color = Type), size = 0.5) +
  labs(x = bquote(Expected -log[10]~italic((P))),
       y = bquote(Observed -log[10]~italic((P))),
       color = NULL) +
  scale_color_manual(
    values = c("#85192d", "#0028a1", "#e89c31"),
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme_my_pub() +
  theme(
    panel.grid.major.x = element_line(linetype = "dotted"),
    panel.grid.major.y = element_line(linetype = "dotted"),
    
    legend.background = element_blank(),
    legend.position = "right",
    
    axis.title = element_text(size = 14), 
    axis.text = element_text(colour = "black", size = 12)
  )

base_plt


# save plt
ggsave(
  "new_Sigma/qq_nullp.pdf",
  base_plt,
  height = 3, width = 4.5,
  useDingbats = TRUE
)

