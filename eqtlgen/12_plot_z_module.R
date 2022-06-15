##############################################
########### plot z-scores of (SNPs, one module) in DGN and eQTLGen ###########
########### for SNPs are eQTLGen specific signals not in DGN ###########
##############################################

**Add statistic summary to check if there is difference between two distributions?**



rm(list = ls())
library(data.table)
library(tidyverse)

source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
snp <- "3:56849749"
chr <- 3


file_eqtlgen_sig_dgn <- 'postanalysis/eqtlgen_spec_to_dgn.txt'
file_snp_meta <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'


# read files -----
eqtlgen_sig_dgn <- fread(file_eqtlgen_sig_dgn, header = TRUE)
snp_meta <- fread(file_snp_meta)
gene_meta <- fread(file_gene_meta)



# determine which modules the given snp correspond to
eqtlgen_sig_dgn <- eqtlgen_sig_dgn %>%
  filter(meta == snp) %>%
  separate(module, into = c(NA, "module"), sep = "module", convert = TRUE)
module_list <- eqtlgen_sig_dgn$module
module_list <- sapply(module_list, function(x) strsplit(x, "module") %>% unlist() %>% .[2] %>% as.numeric()) %>% sort()

file_z_eqtlgen <- paste0("/project2/xuanyao/llw/eQTLGen/z_dgn_166_module/z.module", module_list, ".chr", chr, ".txt.gz")
file_z_dgn <- paste0("/project2/xuanyao/llw/DGN_no_filter_on_mappability/z/z.module", module_list, ".chr", chr, ".txt.gz")



# read z_eqtlgen and z_dgn for specified snp across modules -----
snp_rs <- snp_meta %>% filter(meta == !!snp) %>% pull(SNP)
z_eqtlgen <- lapply(file_z_eqtlgen, function(x){
  tmp_m = basename(x) %>%
    str_split("[.]") %>% unlist() %>% .[2] %>%
    str_split("module") %>% unlist() %>% .[2] %>%
    as.numeric()
  
  fread(x, header = TRUE) %>%
    filter(snp == !!snp_rs) %>%
    pivot_longer(-snp, names_to = "gene", values_to = "z") %>%
    mutate("module" = tmp_m)
  
}) %>%
  rbindlist()

z_dgn <- lapply(file_z_dgn, function(x){
  tmp_m = basename(x) %>%
    str_split("[.]") %>% unlist() %>% .[2] %>%
    str_split("module") %>% unlist() %>% .[2] %>%
    as.numeric()
  
  
  fread(x, header = TRUE) %>%
    filter(snp == !!snp) %>%
    pivot_longer(-snp, names_to = "gene", values_to = "z") %>%
    mutate("module" = tmp_m)
  
}) %>%
  rbindlist()



# convert gene ids -----
gene_meta <- gene_meta %>% separate(Geneid, into = c("gene", NA), sep = "[.]")
z_eqtlgen <- left_join(z_eqtlgen,
                       gene_meta %>% select(gene, GeneSymbol),
                       by = c("gene")) %>%
  select(-gene, -snp) %>%
  rename("gene" = "GeneSymbol") %>%
  mutate("snp" = !!snp)

z_dgn <- z_dgn %>% filter(gene %in% z_eqtlgen$GeneSymbol)



# prep plot data -----
# combine z's across two datasets
plt_df <- rbind(
  mutate(z_dgn, "type" = "DGN"),
  
  mutate(z_eqtlgen, "type" = "eQTLGen")
)

# add annotation of if eQTLGen signal is close to a DGN signal or not
plt_df <- plt_df %>%
  left_join(eqtlgen_sig_dgn %>% select(module, meta, starts_with("if_spec")),
            by = c("module", "snp" = "meta"))

# assign group to each module, for facet wrap. Each row has 10 modules
plt_df$g <- cut(plt_df$module, breaks = c(0, module_list[10*1:20], Inf))
plt_df$module <- factor(plt_df$module)



# Violin plot of modules v.s. z-scores of the given SNP in two datasets -----
base_plt <- ggplot(plt_df, aes(x = module, y = z)) +
  facet_wrap(~g, ncol = 1, scales = "free") +
  geom_violin(aes(fill = type)) +
  geom_point(
    aes(color = type),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 1),
    shape = 16, size = 1, alpha = 1
  ) +
  geom_point(
    data = subset(plt_df, !`if_spec_1e+06`),
    aes(x = module, shape = "eQTLGen specific"), y = 0,
    size = 3, color = "#00cc00"
  ) +
  labs(title = paste0(snp, '\n', snp_rs),
       x = "Module", y = "Z-score",
       color = "Dataset", shape = NULL)

base_plt +
  scale_colour_manual(
    breaks = c("DGN", "eQTLGen"),
    values = c("eQTLGen" = "#0028a1", "DGN" = "#85192d"),
    guide = guide_legend(override.aes = list(size = 2, alpha = 1))
  ) +
  scale_fill_manual(
    breaks = c("DGN", "eQTLGen"),
    values = c("eQTLGen" = "#d6dcef", "DGN" = "#e1c7cc"),
    guide = "none"
  ) +
  scale_shape_manual(
    breaks = "eQTLGen specific", values = 18,
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  theme_my_pub() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed", color = "#999999"),
    
    legend.background = element_blank(),
    legend.position = "right",
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    plot.title = element_text(vjust = -1),
    
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(paste0("plot/z_eqtlgen_dgn_snp", snp, ".pdf"), height = 2*nlevels(plt_df$g), width = 10)

