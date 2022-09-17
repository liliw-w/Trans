###############################################################
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


file.p <- list.files(
  '/project2/xuanyao/llw/DGN_w_cell_prop_pc/p',
  "p.module\\d+.chr\\d+.rds", 
  full.names = TRUE
)


p.obs <- rbindlist(lapply(file.p, function(x)
{tmp_y=readRDS(x);
print(x);
ind = tmp_y < 1e-5
if(!is.null(tmp_y) & sum(ind) > 0){
  tmp_y = tmp_y[ind]
  as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)
}
}))

saveRDS(
  p.obs,
  '/project2/xuanyao/llw/DGN_w_cell_prop_pc/postanalysis/p_all_1e5.rds'
)



###############################################################
###############################################################
rm(list = ls())
library(data.table)


signal_p_cutoff <- 0.1/1e+06/117

res <- readRDS('/project2/xuanyao/llw/DGN_w_cell_prop_pc/FDR/p_all_1e5.rds')
colnames(res) <- c('snp', 'p')
res$q <- NA


fwrite(
  filter(res, p < signal_p_cutoff) %>% arrange(p),
  "/project2/xuanyao/llw/DGN_w_cell_prop_pc/FDR/signals.chr.module.bonf.fdr10.txt",
  col.names = FALSE, quote = FALSE, sep = '\t'
)



###############################################################
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('/home/liliw1/Trans/plot/theme_my_pub_hist.R')


signal <- fread(
  "/project2/xuanyao/llw/DGN_w_cell_prop_pc/FDR/signals.chr.module.bonf.fdr10.txt",
  col.names = c('snp', 'p', 'q')
) %>%
  separate(
    snp, c("module", "chr", "pos"), sep = "[:]", remove = FALSE, convert = TRUE
  ) %>%
  separate(
    module, c(NA, "module"), sep = "module", convert = TRUE
  ) %>%
  unite("meta", chr, pos, sep = ":", remove = FALSE)


signal$module <- factor(
  signal$module
)
signal$chr <- factor(
  signal$chr
)

ggplot(signal, aes(x = chr)) +
  geom_bar(stat = "count") +
  labs(x = "Chromosome") +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave(
  "/project2/xuanyao/llw/DGN_w_cell_prop_pc/plot/signal_chr.pdf",
  width = 5, height = 3
)


ggplot(signal, aes(x = module)) +
  geom_bar(stat = "count") +
  labs(x = "Module") +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.text.x = element_text(angle = 45, size = 8, vjust = 1, hjust = 1)
  )

ggsave(
  "/project2/xuanyao/llw/DGN_w_cell_prop_pc/plot/signal_module.pdf",
  width = 5, height = 3
)


module.snp <- signal$snp
snp <- sort(unique(sapply(strsplit(module.snp, ":"), function(x) paste(x[-1], collapse = ":"))))

str(snp)


