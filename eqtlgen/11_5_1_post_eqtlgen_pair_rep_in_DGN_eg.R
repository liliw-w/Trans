##############################################
########### Interesting examples in replication ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
file_rep <- 'postanalysis/rep_eqtlgen_in_dgn.txt'

file_signal_cis_genes <- "postanalysis/signal_rm_infl_ratio_50.txt_cis_genes.txt"

file_plot_p <- "plot/rep_eqtlgen_in_dgn_p.pdf"
file_plot_p_diff <- "plot/rep_eqtlgen_in_dgn_p_diff.pdf"


# read files -----
sig_rep <- fread(file_rep, header = TRUE)
signal_cis_genes <- fread(file_signal_cis_genes, header = TRUE)


# organize data -----
signal_cis_genes$module <- paste0("module", signal_cis_genes$module)


# add cis genes -----
sig_rep <- sig_rep %>%
  left_join(signal_cis_genes %>% select("module", "SNP", starts_with("near")),
            by = c("module", "SNP"))


# signals in eQTLGen but not DGN -----

## 1. those with extremely small eQTLGen p but extremely large DGN p, SNP v.s. p -----
fig_dat <- sig_rep %>%
  select(p.eqtlgen, p.dgn, if_rep) %>%
  arrange(desc(if_rep), p.eqtlgen, p.dgn) %>%
  mutate("x" = row_number()) %>%
  pivot_longer(c(p.eqtlgen, p.dgn), names_to = "p_Type", values_to = "p")
fig_dat[fig_dat$p == 0, "p"] <- 1e-20
fig_dat$if_rep <- factor(fig_dat$if_rep,
                         levels = c(TRUE, FALSE),
                         labels = c("Replicated", "No Replication") )
fig_dat$p_Type <- factor(fig_dat$p_Type,
                         levels = c("p.eqtlgen", "p.dgn"),
                         labels = c("eQTLGen", "DGN") )


ggplot(fig_dat, aes(x = x, y = -log10(p), color = p_Type)) +
  geom_point(size = 0.7) +
  facet_wrap(vars(if_rep), scales = "free_x") +
  labs(x = "SNP", y = bquote(-Log[10]~italic((P))), color = "Dataset") +
  
  scale_colour_manual(
    values = c("eQTLGen" = "#85192d", "DGN" = "#dd9933"),
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size = unit(0.5, "cm"),
    
    axis.line = element_line(colour = "black"),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text = element_text(colour = "black", size = 14),
    axis.title.y = element_text(vjust = 2, size = 16),
    axis.title.x = element_text(vjust = -0.2, size = 16),
    
    plot.margin = unit(c(10,5,5,5),"mm")
  )
ggsave(file_plot_p, width = 6, height = 4)


## 2. those with extremely small eQTLGen p but extremely large DGN p, difference in p across modules -----
fig_dat <- sig_rep %>%
  mutate("p_diff" = p.dgn - p.eqtlgen) %>%
  select(p_diff, module, if_rep) %>%
  separate(module, c(NA, "module"), "module", convert = TRUE)
fig_dat$module <- factor(fig_dat$module, 1:1000)
fig_dat$if_rep <- factor(fig_dat$if_rep,
                         levels = c(TRUE, FALSE),
                         labels = c("Replicated", "No Replication") )


ggplot(fig_dat, aes(x = module, y = log10(p_diff))) +
  geom_jitter(size = 0.7, color = "#ccac00") +
  facet_wrap(vars(if_rep), scales = "free_x") +
  labs(x = "Module", y = bquote(Log[10] ~Delta ~italic( P)) ) +
  
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size= unit(0.5, "cm"),
    
    axis.line = element_line(colour = "black"),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 6, angle = 90),
    axis.title.y = element_text(vjust = 2, size = 16),
    axis.title.x = element_text(vjust = -0.2, size = 16),
    
    plot.margin = unit(c(10,5,5,5),"mm")
  )
ggsave(file_plot_p_diff, width = 6, height = 3.5)



## 3. those near an interesting gene  -----
sig_rep %>% filter(p.eqtlgen == 0 & p.dgn == 1) %>% View()

sig_rep %>% filter(p.eqtlgen == 0 & p.dgn == 1 & meta == "5:131813219") %>% View()

tmp1 <- sig_rep %>% filter(meta == "5:131813219")
View(tmp1)

tmp2 <- sig_rep %>% filter(SNPChr == 7)
View(tmp2)


sig_rep %>% filter(SNPChr == 5 & grepl("^131\\d{6}", SNPPos) ) %>% View()
sig_rep %>% filter(SNPChr == 5 & grepl("^131\\d{6}", SNPPos) ) %>% count(meta) %>% View()
sig_rep %>% filter(SNPChr == 5 & grepl("^131\\d{6}", SNPPos) ) %>% count(module) %>% View()


sig_rep %>% filter(SNPChr == 3) %>% distinct(meta, .keep_all = TRUE) %>% View()

file_dgn_sig_cis <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/signal_cis_genes.txt'
dgn_sig_cis <- fread(file_dgn_sig_cis)


tmp3 <- dgn_sig_cis %>% filter(nearest_gene == "ARHGEF3")
View(tmp3)

tmp4 <- dgn_sig_cis %>% filter(chr == 7)
View(tmp4)

sig_rep %>%
  filter(SNPChr == 7 & grepl("^50\\d{6}$", SNPPos) ) %>%
  select(module, SNPChr, SNPPos, meta, p.eqtlgen, p.dgn, if_rep) %>%
  separate(module, c(NA, "module"), "module", convert = TRUE) %>%
  count(SNPPos)



### e.g. 1 - loci near gene *IKZF1* -----
fig_dat <- sig_rep %>%
  filter(SNPChr == 7 & grepl("^50\\d{6}$", SNPPos) ) %>%
  select(module, SNPChr, SNPPos, meta, p.eqtlgen, p.dgn, if_rep) %>%
  separate(module, c(NA, "module"), "module", convert = TRUE) %>%
  pivot_longer(c(p.eqtlgen, p.dgn), names_to = "p_Type", values_to = "p")
fig_dat[fig_dat$p == 0, "p"] <- 1e-20

fig_dat$module <- fct_inseq(factor(fig_dat$module))
fig_dat$meta <- fct_infreq(factor(fig_dat$meta))

ggplot(fig_dat, aes(x = module, y = meta, size = -log10(p), color = p_Type)) +
  geom_point(alpha = 1) +
  geom_point(data = filter(fig_dat, if_rep),
             aes(shape = 8),
             size = 3, color = "#0000ff") +
  labs(x = "Module", y = "SNP", size = bquote(~bold( -Log[10] ~italic((P)) ) ), color = "Dataset") +
  
  scale_shape_identity(label = "Replication", guide = guide_legend()) +
  
  scale_colour_manual(
    breaks = c("p.eqtlgen", "p.dgn"),
    labels = c("eQTLGen", "DGN"),
    values = c("p.eqtlgen" = "#771628", "p.dgn" = "#eecc99")
  ) +
  
  scale_size(range = c(0, 5)) +
  
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = NA, fill = NA, size = 1),
    panel.grid.major = element_line(linetype = "dashed", size = 0.3),
    
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    #legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size = unit(0.5, "cm"),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text = element_text(colour = "black", size = 9),
    axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5),
    axis.title.y = element_text(vjust = 2, size = 16),
    axis.title.x = element_text(vjust = -0.2, size = 16),
    
    plot.margin = unit(c(10,5,5,5),"mm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 3)),
         size = guide_legend(order = 2),
         shape = guide_legend(order = 3, title = NULL))
ggsave("plot/eg_snp_near_IKZF1.pdf", width = 11, height = 4)




### e.g. 2 - loci near gene *IRF1* -----
fig_dat <- sig_rep %>%
  filter(SNPChr == 5 & grepl("^131\\d{6}$", SNPPos) ) %>%
  select(module, SNPChr, SNPPos, meta, p.eqtlgen, p.dgn, if_rep) %>%
  separate(module, c(NA, "module"), "module", convert = TRUE) %>%
  pivot_longer(c(p.eqtlgen, p.dgn), names_to = "p_Type", values_to = "p")
fig_dat[fig_dat$p == 0, "p"] <- 1e-20

fig_dat$module <- fct_infreq(factor(fig_dat$module))


ggplot(fig_dat, aes(x = module, y = meta, size = -log10(p), color = p_Type)) +
  geom_point(alpha = 1) +
  labs(x = "Module", y = "SNP", size = bquote(~bold( -Log[10] ~italic((P)) ) ), color = "Dataset") +
  
  scale_colour_manual(
    breaks = c("p.eqtlgen", "p.dgn"),
    labels = c("eQTLGen", "DGN"),
    values = c("p.eqtlgen" = "#771628", "p.dgn" = "#eecc99"),
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  scale_size(range = c(0, 6)) +
  
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = NA, fill = NA, size = 1),
    panel.grid.major = element_line(linetype = "dashed", size = 0.3),
    
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    #legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size = unit(0.5, "cm"),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text = element_text(colour = "black", size = 10),
    axis.text.x = element_text(colour = "black", size = 12, angle = 90, vjust = 0.5),
    axis.title.y = element_text(vjust = 2, size = 16),
    axis.title.x = element_text(vjust = -0.2, size = 16),
    
    plot.margin = unit(c(10,5,5,5),"mm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
ggsave("plot/eg_snp_near_IRF1.pdf", width = 5, height = 6)




### e.g. 3 - Two SNPs `3:56849749` and `3:56865776` near gene *ARHGEF3* -----
fig_dat <- sig_rep %>%
  filter(meta %in% c("3:56849749", "3:56865776") ) %>%
  select(module, meta, p.eqtlgen, p.dgn, if_rep) %>%
  separate(module, c(NA, "module"), "module", convert = TRUE) %>%
  pivot_longer(c(p.eqtlgen, p.dgn), names_to = "p_Type", values_to = "p")
fig_dat$module <- factor(fig_dat$module, 1:1000)
fig_dat[fig_dat$p == 0, "p"] <- 1e-20


ggplot(fig_dat, aes(x = module, y = -log10(p), color = p_Type) ) +
  geom_point(size = 0.7) +
  facet_wrap(vars(meta), scales = "free_x", ncol = 1) +
  labs(x = "Module", y = bquote(-Log[10]~italic((P))), color = "Dataset", title = bquote("Loci near"~italic(ARHGEF3)) ) +
  scale_colour_manual(
    breaks = c("p.eqtlgen", "p.dgn"),
    labels = c("eQTLGen", "DGN"),
    values = c("p.eqtlgen" = "#85192d", "p.dgn" = "#dd9933"),
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size = unit(0.5, "cm"),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", size = 8, angle = 90),
    axis.title.y = element_text(vjust = 2, size = 16),
    axis.title.x = element_text(vjust = -0.2, size = 16),
    
    plot.margin = unit(c(10,5,5,5),"mm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
ggsave("plot/eg_snp_near_ARHGEF3.pdf", width = 6, height = 5)



### p values in DGN of signals for module4 in loci ARHGEF3
fig_dat <- dgn_sig_cis %>%
  filter(chr == 3 & module == "module4") %>%
  select(pos, p) %>%
  arrange(pos) %>%
  mutate("M" = "Module 4")
fig_dat[fig_dat$p == 0, "p"] <- 1e-20

ggplot(fig_dat, aes(x = factor(pos), y = -log10(p)) ) +
  geom_point(size = 1.5, color = "#dd9933") +
  geom_vline(xintercept = "56849749", linetype = "dashed", color = "grey") +
  geom_vline(xintercept = "56865776", linetype = "dashed", color = "grey") +
  facet_wrap(vars(M)) +
  labs(x = "Chromosome 3", y = bquote(-Log[10]~italic((P))) ) +
  
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size = unit(0.5, "cm"),
    
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", size = 8, angle = 90),
    axis.title.y = element_text(vjust = 2, size = 16),
    axis.title.x = element_text(vjust = -0.2, size = 16),
    
    plot.margin = unit(c(10,5,5,5),"mm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
ggsave("plot/eg_snp_near_ARHGEF3_module4_DGN.pdf", width = 6, height = 3.5)


