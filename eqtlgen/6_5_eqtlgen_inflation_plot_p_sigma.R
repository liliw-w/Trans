###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Compare p's using real-data sigma and nullz sigma ###########
###############################################################
module <- as.numeric(snakemake@params[['module']])

# load packages -----
library(tidyverse)
source('~/Trans/plot/plt_qq_multi.R')


# I/O & paras -----
thre_p_z <- 1e-4

file_Sigma <- paste0("inflation/Sigma_DGN_module", module, ".rds")
file_p <- paste0('inflation/p_', basename(file_Sigma))
file_p_nullz_seq <- list.files(path = "inflation",
                               pattern = paste0("^p_Sigma_nullz\\d+\\_", basename(file_Sigma)),
                               full.names = TRUE)
# for formatting plt tile
file_coexp_module <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
file_null_SNP <- 'null_SNP/num_nullSNP.rds'

## output -----
file_plot_hist <- paste0('inflation/plot/plot_hist_', basename(file_Sigma), '.pdf')
file_plot_qq <- paste0('inflation/plot/plot_qq_', basename(file_Sigma), '.pdf')


# read files & important paras -----
## for formatting plt tile
coexp_module <- readRDS(file_coexp_module)$moduleLabels
res_nullSNP <- readRDS(file_null_SNP)

## all p from various sigma
p_all <- lapply(c(file_p, file_p_nullz_seq), readRDS) %>% do.call(cbind, .) %>% as_tibble()


# organize data -----
## extract numbers from file names to calculate actual ratio used -----
digit_in_file <- strsplit(file_p_nullz_seq, "\\D") %>%
  unlist() %>%
  as.numeric() %>%
  .[!is.na(.)] %>%
  matrix(ncol = length(file_p_nullz_seq)) %>%
  as_tibble()
n_nullz_seq <- digit_in_file[1, ] %>% unlist(use.names = FALSE)
K <- sum(coexp_module==module) #digit_in_file[3, ] %>% unlist(use.names = FALSE) %>% unique()

ratio <- n_nullz_seq / K


## assign names to p's with their correspoding sigma -----
#colnames(p_all)[-1] <- paste0("Nnull", n_nullz_seq, "_Ratio", ratio)
colnames(p_all)[-1] <- paste0("Ratio=", ratio)
colnames(p_all)[1] <- paste0("Real Sigma")

## format plt title -----
eqtlgen_ratio <- res_nullSNP %>%
  filter(thre_z == thre_p_z) %>%
  mutate(prop_dim = num_nullSNP_indep/module_size) %>%
  select(num_nullSNP_indep, module, module_size, prop_dim) %>%
  filter(module == {{module}}) %>%
  pull(prop_dim) %>%
  round(digits = 2)
plt_title <- paste0("M", module, "_Size", K, "_Ratio", eqtlgen_ratio)



# visualization -----
## 1. hist -----
plt_dt <- p_all %>% pivot_longer(everything(), names_to = "Sigma", values_to = "p")
plt_dt$Sigma <- factor(plt_dt$Sigma,
                       levels = colnames(p_all)[c(1, order(ratio, decreasing = TRUE) + 1)])

plt_obj <- ggplot(plt_dt, aes(x = p)) +
  geom_histogram(binwidth = 0.01, fill = "#3a3a3a") +
  facet_wrap(vars(Sigma), scales = "free_y") +
  labs(x = "P-value", y = "Number of SNPs", title = plt_title) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.6, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        
        plot.title = element_text(hjust = 0.5, face = "bold"),
        
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "#f5f5f5", colour = "black", size = 0.6, linetype = "solid")
  )
ggsave(file_plot_hist, plt_obj, width = 9, height = 6)
saveRDS(plt_obj, paste0(file_plot_hist, ".rds"))


## 2. QQ-plot -----
fig <- qqplot(p_all,
              is_group_numerical_order = TRUE,
              group_title = "Sigma",
              group_order = colnames(p_all)[c(1, order(ratio, decreasing = TRUE) + 1)])
plt_obj <- fig +
  labs(title = plt_title) +
  theme(text = element_text(size = 6),
        plot.title = element_text(size = 16),
        legend.position = "none",
        legend.title = element_text(angle = 90),
        legend.text = element_text(size = 10)) +
  guides(shape = guide_legend(ncol = 3, override.aes = list(size = 3)))

ggsave("tmp.pdf", plt_obj, width = 4, height = 4, useDingbats = TRUE)


ggsave(file_plot_qq, plt_obj, width = 6, height = 6, useDingbats = TRUE)
saveRDS(plt_obj, paste0(file_plot_qq, ".rds"))


