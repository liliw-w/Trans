##############################################
############ archie Empirical p and power ############
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')

theme_set(
  theme_my_pub(legend.position = "right")
)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    archie_cc_value_z.rds
    archie_selected_gene_z.rds
    archie_null_cc_value_z.rds
    archie_null_selected_gene_z.rds
    0.95
    0.05
    archie_p_sim.pdf
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_cc_alt <- args[1]
file_selected_gene_alt <- args[2]
file_cc_null <- args[3]
file_selected_gene_null <- args[4]
ci_level <- as.numeric(args[5])
fdr_level <- as.numeric(args[6])

## output -----
file_fig <- args[7]


# read files -----
cc_alt <- readRDS(file_cc_alt)
selected_gene_alt <- readRDS(file_selected_gene_alt)
cc_null <- readRDS(file_cc_null)
selected_gene_null <- readRDS(file_selected_gene_null)

cc_alt <- cc_alt[,,,, drop = TRUE]
cc_null <- cc_null[,, drop = TRUE]
cc_null <- cc_null[!is.na(cc_null)]


# empirical p -----
n_null <- length(cc_null)

emp_p <- apply(
  cc_alt,
  c(1, 2),
  function(x) sum(cc_null > x)
) / n_null

# emp_p <- apply(
#   cc_alt,
#   1,
#   function(x) qvalue::empPvals(x, cc_null)
# ) %>%
#   t()


# power -----
## adjust p
adj_p <- apply(
  emp_p,
  1,
  function(x) qvalue::qvalue(x, fdr.level = fdr_level)$qvalues
) %>%
  t()

## power
adj_power <- apply(adj_p, 1, function(x) sum(x < fdr_level))




# qq plot & histogram & power -----
n_test <- ncol(cc_alt)

# order statistic of null p
expected <- seq(1, n_test) / (n_test+1)
lexp <- -log10(expected)

ci_l <- -log10( qbeta(p = (1 - ci_level) / 2, shape1 = 1:n_test, shape2 = n_test:1) )
ci_r <- -log10( qbeta(p = (1 + ci_level) / 2, shape1 = 1:n_test, shape2 = n_test:1) )

# obs
observed <- apply(t(emp_p), 2, sort) %>% as.data.frame()
lobs <- -log10(observed)


# data for plt
df_plt = cbind(data.frame(x = lexp, ci_l = ci_l, ci_r = ci_r), lobs) %>%
  pivot_longer(-c(x, ci_l, ci_r), names_to = "Type", values_to = "y")


ggpubr::ggarrange(
  # 1. qq fig
  ggplot(df_plt, aes(x = x, y = y, group = Type)) +
    geom_ribbon(aes(ymin = ci_l, ymax = ci_r), fill = "#e5e5e5", color = "#e5e5e5") +
    geom_abline(slope = 1, intercept = 0, color = "#595959", size = 0.7) +
    geom_point(aes(color = Type), size = 0.5) +
    labs(
      x = bquote(Expected -log[10]~italic((P))),
      y = bquote(Observed -log[10]~italic((P)))
    ) +
    scale_color_manual(
      values = c(
        "#85192d", "#0028a1", "#e89c31",
        RColorBrewer::brewer.pal(8, "Dark2"),
        RColorBrewer::brewer.pal(8, "Set1")
      ) 
    ),
  
  
  # 2. hist fig of p
  as_tibble(emp_p, rownames = "model") %>%
    pivot_longer(
      cols = !model,
      names_to = NULL, values_to = "p"
    ) %>%
    ggplot() +
    facet_wrap(~model) +
    geom_histogram(
      aes(x = p)
    ) +
    labs(x = "Empirical P", y = "Simulation"),
  
  
  # 3. power
  as_tibble(adj_p, rownames = "model") %>%
    pivot_longer(
      cols = !model,
      names_to = NULL, values_to = "p"
    ) %>%
    ggplot() +
    geom_boxplot(aes(x = model, y = -log10(p))) +
    geom_hline(yintercept = -log10(fdr_level), linetype = "dashed", color = "maroon") +
    labs(y = bquote(-log[10]~italic((Adj_P)))),
  
  nrow = 3, labels = LETTERS[1:3], heights = c(1, 1.4, 1)
)



# print out key message or write out -----
ggsave(
  file_fig,
  height = 6, width = 4
)

