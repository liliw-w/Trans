##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggrepel)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
ratio <- 50
#dis_cis <- 1e+5


#file_qtl <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')
file_signal_cis_genes <- paste0('postanalysis/signal_cis_genes_rm_infl_ratio_', ratio, '.txt')
file_TF <- '/scratch/midway2/liliw1/figures/Manhattan/TF_names_v_1.01.txt'
file_chr_pos <- '/scratch/midway2/liliw1/sig_module_chr/chromosome_location.rds'

## output -----
file_manhattan_gene


# read files -----
signal_cis_genes <- fread(file_signal_cis_genes, header = TRUE)
TF <- fread(file_TF, header = FALSE, col.names = "TF_gene")
chr_pos <- readRDS(file_chr_pos)


# how many signal pairs have TF as cis- genes -----
tmp <- sapply(TF$TF_gene, function(x) grep(x, signal_cis_genes$near_genes) )
tmp_df <- tibble(TF_gene = names(tmp),
                 n_signal = sapply(tmp, length) ) %>%
  arrange(desc(n_signal))
View(tmp_df)



# 0. one p for one snp across modules
snp_cis_genes <- signal_cis_genes %>%
  group_by(meta, SNPChr, SNPPos, nearest_gene, near_genes, nearest_dis, near_dis) %>%
  summarise('p' = min(p),
            'n_module' = n()) %>%
  mutate('-logp' = -log10(p)) %>%
  rename("chr" = "SNPChr", "pos" = "SNPPos") %>%
  ungroup()

snp_cis_genes <- separate_rows(
  snp_cis_genes,
  near_genes, near_dis,
  sep = ";", convert = TRUE
)
  #filter(
  #  nearest_dis < dis_cis/2 & near_dis < dis_cis/2
  #)

filter(
  snp_cis_genes,
  nearest_gene %in% !!TF$TF_gene
) %>%
  distinct(nearest_gene) %>%
  View()

# change 0 p -----
snp_cis_genes[is.infinite(snp_cis_genes$`-logp`), "-logp"] <-
  (max(snp_cis_genes$`-logp`[!is.infinite(snp_cis_genes$`-logp`)]) + 1) %>% ceiling()


# 0. gene of interest -----
gene_of_interest <- c("NFKBIA", "PLAGL1", "NFE2", "IKZF1", "KLF1", "KLF14", "NFKB1", "ZNF229", "BAZ2B", #TF
                      "ARHGEF3", "SENP7")


# 1. nearest & near gene of interest col -----
snp_cis_genes <- mutate(
  snp_cis_genes,
  "if_gene_of_interest_nearest" = nearest_gene %in% !!gene_of_interest,
  "if_gene_of_interest_near" = near_genes %in% !!gene_of_interest,
)

plt_dat <- snp_cis_genes %>%
  left_join(chr_pos, by = c("chr" = "CHR")) %>%
  mutate(def_pos = pos + tot)


# 4. same annotation, pick the smallest p ----
annot_gene <- filter(plt_dat, if_gene_of_interest_nearest) %>%
  mutate("label_gene_of_interest" = nearest_gene) %>%
  group_by(label_gene_of_interest) %>%
  summarise(`-logp` = max(`-logp`),
            chr = chr, pos = pos, def_pos = def_pos) %>%
  ungroup() %>%
  distinct(label_gene_of_interest, .keep_all = TRUE)



annot_gene <- filter(plt_dat, if_gene_of_interest_near) %>%
  mutate("label_gene_of_interest" = near_genes) %>%
  group_by(label_gene_of_interest) %>%
  summarise(`-logp` = max(`-logp`),
            chr = chr, pos = pos, def_pos = def_pos) %>%
  ungroup() %>%
  distinct(label_gene_of_interest, .keep_all = TRUE)


annot_snp <- filter(plt_dat, if_gene_of_interest_near) %>%
  mutate("label_gene_of_interest" = near_genes)



# plot -logp for (cht pos, module) -----
ggplot(plt_dat, aes(x = def_pos, y = `-logp`)) +
  geom_rect(
    aes(xmin = tot, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = factor(chr))
  ) +
  geom_point(aes(color = factor(chr)), alpha = 0.5, size = 1, shape = 16) +
  #geom_point(data = subset(don, is_annotate_TF=="yes"),
  #           color = "#b20000", fill = "#b20000", shape=23, size = 3) +
  labs(x = "Chromosome", y = quote(-Log[10](P))) +
  scale_x_continuous(
    limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
    label = chr_pos$CHR,
    breaks = chr_pos$center,
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(6.5, 21) ) +
  scale_color_manual(values = rep(c("#7e7e7e", "#bfbfbf"), 22), guide = "none") +
  scale_fill_manual(values = rep(c("#efefef", "#ffffff"), 22), guide = "none") +
  theme_my_pub() +
  theme(
    #axis.text.x = element_text(angle = 90, size = 8)
  ) +
  geom_point(data = annot_snp,
             color = "#b20000", alpha = 0.5, size = 1) +
  geom_point(data = annot_gene,
             color = "#b20000", fill = "#b20000", shape=23, size = 3) +
  geom_text_repel(data = annot_gene,
                   aes(label = label_gene_of_interest),
                   segment.colour="black",
                   size = 4,
                   min.segment.length = 0,
                   max.overlaps = 5,
                   nudge_x = -0.5,
                   nudge_y = 30,
                   box.padding = 1,
                   segment.curvature = -0.1,
                   segment.ncp = 5,
                   segment.angle = 20,
                   direction = "y",
                   hjust = "left",
                   segment.linetype = 6,
                   arrow = arrow(length = unit(0.015, "npc"))
  )

ggsave('tmp_manhattan_cis_gene.pdf', width = 7, height = 4)

