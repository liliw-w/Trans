###############################################################
########### plot unique significant SNPs ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source("~/Trans/plot/theme_my_pub.R")


file_signal <- 'postanalysis/signal.txt'
n_SNP_total <- 9918


### read data
signal <- fread(file_signal)


### how many signals each module has?
signal %>% count(module) %>% arrange(desc(n))


### how many modules each signal correspond to?
signal %>%
  group_by(SNP, meta, SNPChr, SNPPos) %>%
  summarise("n_module" = n(), "minp" = min(p) ) %>%
  arrange(minp, desc(n_module), SNPChr)


### 1. plot how many signals each module has?
signal$module <- factor(signal$module, levels = 1:1000)

fig <- ggplot(signal, aes(x = module)) +
  geom_bar(stat = "count") +
  geom_hline(yintercept = n_SNP_total,
             color = "#990000", linetype = "dashed", size = 1) +
  labs(x = "Module", y = "Number of signals\n(SNPs)")
fig +
  theme_my_pub() +
  theme(axis.text.x = element_text(size = 6, angle = 90),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"))

ggsave("plot/signal_module.pdf", width = 10, height = 4)


### 2. plt how many signals are on each chr?
signal$SNPChr <- factor(signal$SNPChr, levels = 1:1000)

fig <- ggplot(signal %>% distinct(SNP, .keep_all = TRUE),
              aes(x = SNPChr)) +
  geom_bar(stat = "count") +
  labs(x = "Chromosome", y = "Number of signals\n(SNPs)")
fig +
  theme_my_pub() +
  theme(axis.text.x = element_text(angle = 60),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"))

ggsave("plot/signal_chr.pdf", width = 6, height = 4)


### 3. plt how many signals for each chr and module?
# set module & chr order on the figure
plt_df <- signal %>%
  group_by(module, SNPChr) %>%
  summarise(n_M_chr = n()) %>%
  ungroup() %>%
  arrange(desc(n_M_chr))
plt_df$module = factor(plt_df$module, levels = plt_df$module, labels = plt_df$module)
plt_df$SNPChr = factor(plt_df$SNPChr, levels = (plt_df$SNPChr) |> unique() |> rev())


fig <- ggplot(plt_df, aes(x = module, y = SNPChr)) +
  geom_tile(aes(fill = n_M_chr)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[3:8]) +
  labs(x = "Module", y = "Chromosome", fill = "SNPs")
fig +
  theme_my_pub(legend.text.size = 10, legend.position = "right") +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.ticks.x = element_blank())

ggsave("plot/signal_module_chr.pdf", width = 12, height = 5)
