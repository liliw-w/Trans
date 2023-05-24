##############################################
########### check how many eqtlgen univariate signals are also signals for coexpression modules or pathways ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
file_sig_eqtlgen <- '/project2/xuanyao/llw/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz'
file_sig_coexp <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/signal_uniq_rm_infl_ratio_50.txt'
file_sig_msig <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/signal_uniq_rm_infl_ratio_50.txt'

## output -----
file_plt <- "/project2/xuanyao/llw/eQTLGen/eqtlgen_sig_in_coexp_msig.pdf"


# read files -----
sig_eqtlgen <- fread(file_sig_eqtlgen)
sig_coexp <- fread(file_sig_coexp, header = TRUE)
sig_msig <- fread(file_sig_msig, header = TRUE)


# check if eqtlgen univariate signals are also signals for coexpression modules or pathways -----
sig_eqtlgen_rep <- sig_eqtlgen %>%
  unite(col = "snp_id", SNPChr, SNPPos, sep = ":", remove = FALSE) %>%
  group_by(snp_id) %>%
  summarise(n = n()) %>%
  ungroup()

sig_eqtlgen_rep$group <- cut(sig_eqtlgen_rep$n, c(1:11, Inf), c(1:10, ">10"), right = FALSE)
sig_eqtlgen_rep$if_coexp <- sig_eqtlgen_rep$snp_id %in% sig_coexp$meta
sig_eqtlgen_rep$if_msig <- sig_eqtlgen_rep$snp_id %in% sig_msig$meta
sig_eqtlgen_rep$if_coexp_or_msig <- sig_eqtlgen_rep$if_coexp | sig_eqtlgen_rep$if_msig


# add a line representing transPCO specific signals -----
transpco_trans <- rbind(sig_coexp, sig_msig) %>% distinct(meta)
n_transpco_spec <- nrow(transpco_trans) - sum(transpco_trans$meta %in% sig_eqtlgen_rep$snp_id)

cat(
  'There are', nrow(sig_eqtlgen_rep), "eQTLGen trans, ",
  sum(sig_eqtlgen_rep$if_coexp_or_msig), "are either coexp or msig trans. \n\n",
  
  'There are', nrow(transpco_trans), "coexp or msig trans, ",
  n_transpco_spec, "are not eQTLGen trans. \n\n"
)


# plot if eqtlgen univariate signals are also signals for coexpression modules or pathways -----
plt_dat <- group_by(sig_eqtlgen_rep, group) %>%
  summarise(
    eqtlgen = n(),
    if_msig = sum(if_msig),
    if_coexp = sum(if_coexp),
    if_coexp_or_msig = sum(if_coexp_or_msig)
  ) %>%
  pivot_longer(!group, names_to = "type", values_to = "n_sig")


ggplot(plt_dat, aes(x = group, y = n_sig, fill = type)) +
  #geom_col(
  #  data = subset(plt_dat, type %in% c("eqtlgen", "if_msig", "if_coexp")),
  #  position = position_dodge()
  #) +
  geom_col(
    data = subset(plt_dat, type %in% c("eqtlgen", "if_coexp_or_msig")),
    position = "identity"
    #position = position_nudge(x = -0.3),
    #width = 0.3
  ) +
  geom_hline(yintercept = n_transpco_spec, linetype = "dashed") +
  labs(x = "Number of trans target genes", y = "Number of trans-eQTLs", fill = "Signals") +
  scale_fill_manual(
    breaks = c("eqtlgen", "if_coexp_or_msig"),
    values = c(
      "eqtlgen" = "#A6CEE3", 
      "if_coexp_or_msig" = "#1F78B4", 
      "if_coexp" = "#743535", 
      "if_msig" = "#cf8825"
    ),
    labels = c(
      "eqtlgen" = "eQTLGen", 
      "if_coexp_or_msig" = "Coexpression module or Pathway"
    )
  ) +
  theme_my_pub(
    base_size = 12,
    legend.text.size = 8,
    axis.title.size = 12,
    axis.text.size = 10
  ) +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    legend.background = element_blank(),
    legend.position = "right",
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(file_plt, width = 5.5, height = 3)

