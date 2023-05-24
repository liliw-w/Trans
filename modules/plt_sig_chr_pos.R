##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
file_chr_pos <- '/scratch/midway2/liliw1/sig_module_chr/chromosome_location.rds'
file_qtl <- '/project2/xuanyao/llw/MODULES/MSigDB/FDR/signals.chr.module.perm10.fdr10.txt'
file_coexp_module <- '/project2/xuanyao/llw/MODULES/MSigDB/result/coexp.module.rds'


## output -----
file_plt_signal_module_chr_pos <- '/project2/xuanyao/llw/MODULES/MSigDB/plot/signal_module_chr_chr_pos.pdf'


# read files -----
chr_pos <- readRDS(file_chr_pos)
qtl <- fread(file_qtl, header = FALSE, col.names = c("signal", "p", "q"))
coexp_module <- readRDS(file_coexp_module)


# organize data -----
qtl <- qtl %>%
  separate(signal, c("module", "chr", "pos"), ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_ID", chr, pos, sep = ":", remove = FALSE) %>%
  separate(module, c(NA, "module"), "module", convert = TRUE)


#qtl <- qtl %>%
#  rename("chr" = "SNPChr", "pos" = "SNPPos")

# add MSigDB pathway as module annotation -----
coexp_module_annot <- enframe(
  coexp_module$moduleName,
  name = "annot_module",
  value = "module"
)


qtl <- left_join(
  qtl, coexp_module_annot,
  by = "module"
)


plt_dat <- qtl %>%
  left_join(chr_pos, by = c("chr" = "CHR")) %>%
  mutate(def_pos = pos + tot) %>%
  arrange(module)
plt_dat$module <- factor(plt_dat$module, levels = plt_dat$module, labels = plt_dat$annot_module)


# plot -logp for (cht pos, module) -----
ggplot(plt_dat) +
  geom_rect(
    data = chr_pos,
    aes(xmin = tot, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = factor(CHR))
  ) +
  geom_hline(
    aes(yintercept = module),
    linetype = "dotted", color = "#e5e5e5") +
  geom_point(
    aes(
      x = def_pos, y = module,
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
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.ticks.y = element_blank()
  )


# print out key message or write out -----
ggsave(
  file_plt_signal_module_chr_pos,
  width = 8, height = 6
)

