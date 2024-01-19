##############################################
########### Plot to check if drug targets are more likely to be trans signals of gene modules ###########
########### or if trans signals of gene modules are more likely to be drug targets ###########
##############################################
rm(list = ls())
library(tidyverse)


# 0. prep data for plotting -----
## input
file_target_enrich <- list.files(".", "*target_enrich_hallmark.txt", full.names = TRUE)

## output
file_fig_target_enrich_target_vs_trans_snp <- "target_enrich_target_vs_trans_snp_hallmark.pdf"
file_fig_target_enrich_target_vs_trans_gene <- "target_enrich_target_vs_trans_gene_hallmark.pdf"



## read data -----
target_enrich <- lapply(
  file_target_enrich,
  FUN = data.table::fread
) %>%
  bind_rows()



## use snps for plotting -----
plt_dat <- select(
  target_enrich,
  trait, if_trans_sig, if_target_all, if_target_trait, if_immune_blood_anno
) %>%
  count(trait, if_trans_sig, if_target_all, if_target_trait, if_immune_blood_anno)
plt_dat$if_target_all <- factor(plt_dat$if_target_all, levels = c('FALSE', 'TRUE'))
plt_dat$if_target_trait <- factor(plt_dat$if_target_trait, levels = c('FALSE', 'TRUE'))

plt_dat$if_trans_sig_immune_blood_anno <- plt_dat$if_trans_sig & as.logical(plt_dat$if_immune_blood_anno)
plt_dat$if_trans_sig <- factor(plt_dat$if_trans_sig, levels = c('TRUE', 'FALSE'))
plt_dat$if_trans_sig_immune_blood_anno <- factor(plt_dat$if_trans_sig_immune_blood_anno, levels = c('TRUE', 'FALSE'))


## use genes for plotting -----
plt_dat_gene <- select(
  target_enrich,
  trait, nearest_gene, if_trans_sig, if_target_all, if_target_trait, if_immune_blood_anno
) %>%
  group_by(trait, nearest_gene) %>%
  summarise(
    if_trans_sig = any(if_trans_sig),
    if_target_all = any(if_target_all),
    if_target_trait = any(if_target_trait),
    if_immune_blood_anno = any(if_immune_blood_anno)
  ) %>%
  ungroup() %>%
  count(trait, if_trans_sig, if_target_all, if_target_trait, if_immune_blood_anno)

plt_dat_gene$if_target_all <- factor(plt_dat_gene$if_target_all, levels = c('FALSE', 'TRUE'))
plt_dat_gene$if_target_trait <- factor(plt_dat_gene$if_target_trait, levels = c('FALSE', 'TRUE'))

plt_dat_gene$if_trans_sig_immune_blood_anno <- plt_dat_gene$if_trans_sig & as.logical(plt_dat_gene$if_immune_blood_anno)
plt_dat_gene$if_trans_sig <- factor(plt_dat_gene$if_trans_sig, levels = c('TRUE', 'FALSE'))
plt_dat_gene$if_trans_sig_immune_blood_anno <- factor(plt_dat_gene$if_trans_sig_immune_blood_anno, levels = c('TRUE', 'FALSE'))




# # 1. trans v.s. target -----
# 
# ## 1.1. use snps -----
# ggpubr::ggarrange(
#   ggplot(plt_dat) +
#     facet_wrap(~trait, nrow = 1) +
#     geom_col(
#       aes(
#         x = if_trans_sig,
#         y = n,
#         fill = if_target_all
#       )
#     ) +
#     labs(
#       x = "If trans signal",
#       y = "Number of signals",
#       fill = "If target \n(all traits)"
#     ),
#   
#   ggplot(plt_dat) +
#     facet_wrap(~trait, nrow = 1) +
#     geom_col(
#       aes(
#         x = if_trans_sig,
#         y = n,
#         fill = if_target_trait
#       )
#     ) +
#     labs(
#       x = "If trans signal",
#       y = "Number of signals",
#       fill = "If target \n(corresponding trait)"
#     ),
#   
#   ggplot(plt_dat) +
#     facet_wrap(~trait, nrow = 1) +
#     geom_col(
#       aes(
#         x = if_trans_sig & as.logical(if_immune_blood_anno),
#         y = n,
#         fill = if_target_all
#       )
#     ) +
#     labs(
#       x = "If trans signal of immune-related modules",
#       y = "Number of signals",
#       fill = "If target \n(all traits)"
#     ),
#   
#   
#   ncol = 1, labels = LETTERS[1:3],
#   align = "hv"
# )
# 
# 
# ggsave(
#   filename = 'target_enrich_trans_vs_target_snp.pdf',
#   width = 8, height = 6
# )





# 2. target v.s. trans -----

## 2.1. use snps -----
ggpubr::ggarrange(
  ggplot(plt_dat) +
    facet_wrap(~trait, nrow = 1) +
    geom_col(
      aes(
        x = if_target_all,
        y = n,
        fill = if_trans_sig
      )
    ) +
    geom_text(
      data = group_by(
        plt_dat,
        trait, if_target_all, if_trans_sig
      ) %>%
        summarise(
          n = sum(n)
        ) %>%
        arrange(
          if_target_all, if_trans_sig = desc(if_trans_sig),
          .by_group = TRUE
        ) %>%
        mutate(label_y = cumsum(n)) %>%
        ungroup(),
      aes(
        x = if_target_all,
        y = label_y,
        label = n
      ),
      size = 2
    ) +
    labs(
      x = "If target (all traits)",
      y = "Number of SNPs",
      fill = "If trans signal"
    ),
  
  ggplot(plt_dat) +
    facet_wrap(~trait, nrow = 1) +
    geom_col(
      aes(
        x = if_target_trait,
        y = n,
        fill = if_trans_sig
      )
    ) +
    geom_text(
      data = group_by(
        plt_dat,
        trait, if_target_trait, if_trans_sig
      ) %>%
        summarise(
          n = sum(n)
        ) %>%
        arrange(
          if_target_trait, if_trans_sig = desc(if_trans_sig),
          .by_group = TRUE
        ) %>%
        mutate(label_y = cumsum(n)) %>%
        ungroup(),
      aes(
        x = if_target_trait,
        y = label_y,
        label = n
      ),
      size = 2
    ) +
    labs(
      x = "If target (corresponding trait)",
      y = "Number of SNPs",
      fill = "If trans signal"
    ),
  
  ggplot(plt_dat) +
    facet_wrap(~trait, nrow = 1) +
    geom_col(
      aes(
        x = if_target_all,
        y = n,
        fill = if_trans_sig_immune_blood_anno
      )
    ) +
    geom_text(
      data = group_by(
        plt_dat,
        trait, if_target_all, if_trans_sig_immune_blood_anno
      ) %>%
        summarise(
          n = sum(n)
        ) %>%
        arrange(
          if_target_all, if_trans_sig_immune_blood_anno = desc(if_trans_sig_immune_blood_anno),
          .by_group = TRUE
        ) %>%
        mutate(label_y = cumsum(n)) %>%
        ungroup(),
      aes(
        x = if_target_all,
        y = label_y,
        label = n
      ),
      size = 2
    ) +
    labs(
      x = "If target (all traits)",
      y = "Number of SNPs",
      fill = "If trans signal \n of immune-related modules"
    ),
  
  
  ncol = 1, labels = LETTERS[1:3],
  align = "hv"
)

ggsave(
  filename = file_fig_target_enrich_target_vs_trans_snp,
  width = 8, height = 6
)



## 2.2. use genes -----
ggpubr::ggarrange(
  ggplot(plt_dat_gene) +
    facet_wrap(~trait, nrow = 1) +
    geom_col(
      aes(
        x = if_target_all,
        y = n,
        fill = if_trans_sig
      )
    ) +
    geom_text(
      data = group_by(
        plt_dat_gene,
        trait, if_target_all, if_trans_sig
      ) %>%
        summarise(
          n = sum(n)
        ) %>%
        arrange(
          if_target_all, if_trans_sig = desc(if_trans_sig),
          .by_group = TRUE
        ) %>%
        mutate(label_y = cumsum(n)) %>%
        ungroup(),
      aes(
        x = if_target_all,
        y = label_y,
        label = n
      ),
      size = 2
    ) +
    labs(
      x = "If target (all traits)",
      y = "Number of nearest genes",
      fill = "If trans signal"
    ),
  
  ggplot(plt_dat_gene) +
    facet_wrap(~trait, nrow = 1) +
    geom_col(
      aes(
        x = if_target_trait,
        y = n,
        fill = if_trans_sig
      )
    ) +
    geom_text(
      data = group_by(
        plt_dat_gene,
        trait, if_target_trait, if_trans_sig
      ) %>%
        summarise(
          n = sum(n)
        ) %>%
        arrange(
          if_target_trait, if_trans_sig = desc(if_trans_sig),
          .by_group = TRUE
        ) %>%
        mutate(label_y = cumsum(n)) %>%
        ungroup(),
      aes(
        x = if_target_trait,
        y = label_y,
        label = n
      ),
      size = 2
    ) +
    labs(
      x = "If target (corresponding trait)",
      y = "Number of nearest genes",
      fill = "If trans signal"
    ),
  
  ggplot(plt_dat_gene) +
    facet_wrap(~trait, nrow = 1) +
    geom_col(
      aes(
        x = if_target_all,
        y = n,
        fill = if_trans_sig_immune_blood_anno
      )
    ) +
    geom_text(
      data = group_by(
        plt_dat_gene,
        trait, if_target_all, if_trans_sig_immune_blood_anno
      ) %>%
        summarise(
          n = sum(n)
        ) %>%
        arrange(
          if_target_all, if_trans_sig_immune_blood_anno = desc(if_trans_sig_immune_blood_anno),
          .by_group = TRUE
        ) %>%
        mutate(label_y = cumsum(n)) %>%
        ungroup(),
      aes(
        x = if_target_all,
        y = label_y,
        label = n
      ),
      size = 2
    ) +
    labs(
      x = "If target (all traits)",
      y = "Number of nearest genes",
      fill = "If trans signal \n of immune-related modules"
    ),
  
  
  ncol = 1, labels = LETTERS[1:3],
  align = "hv"
)

ggsave(
  filename = file_fig_target_enrich_target_vs_trans_gene,
  width = 8, height = 6
)

