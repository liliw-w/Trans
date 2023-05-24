##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
ratio <- 50

file_chr_pos <- '/scratch/midway2/liliw1/sig_module_chr/chromosome_location.rds'
file_qtl <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')

## output -----
file_plt_signal_module_chr_pos <- paste0('plot/signal_module_chr_rm_infl_ratio_', ratio, '_chr_pos.pdf')


# read files -----
chr_pos <- readRDS(file_chr_pos)
qtl <- fread(file_qtl, header = TRUE)

# change signal FDR level to 0.05, to be consistent with the paper -----
thre_p <- 0.05/9918/129
qtl <- filter(qtl, p < thre_p)


# organize data -----
qtl <- qtl %>%
  rename("chr" = "SNPChr", "pos" = "SNPPos")

plt_dat <- qtl %>%
  left_join(chr_pos, by = c("chr" = "CHR")) %>%
  mutate(def_pos = pos + tot)


# plot -logp for (cht pos, module) -----
ggplot(plt_dat) +
  geom_rect(
    data = chr_pos,
    aes(xmin = tot, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = factor(CHR))
  ) +
  geom_hline(
    aes(yintercept = factor(module)),
    linetype = "dotted", color = "#e5e5e5") +
  geom_point(
    aes(
      x = def_pos, y = factor(module),
      size = -log10(p), color = factor(chr)
    ),
    alpha = 0.5, shape = 1) +
  labs(y = "Gene modules", x = "Chromosome", size = quote(-Log[10](P))) +
  scale_x_continuous(
    limits = c(0, max(chr_pos$center)*2 - max(chr_pos$tot)),
    label = chr_pos$CHR,
    breaks = chr_pos$center,
    expand = c(0, 0)
  ) +
  scale_size(guide = guide_legend(override.aes = list(alpha = 1)), range = c(0.5, 3)) +
  scale_color_manual(values = rep(c("#800080", "#daa520"), 22), guide = "none") +
  scale_fill_manual(values = rep(c("#e5e5e5", "#ffffff"), 22), guide = "none") +
  theme_my_pub() +
  theme(
    axis.text.y = element_text(size = 8)
  )


# print out key message or write out -----
ggsave(
  file_plt_signal_module_chr_pos,
  width = 6.5, height = 9.5
)

