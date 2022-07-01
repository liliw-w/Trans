rm(list = ls())
setwd('/scratch/midway2/liliw1/Manhattan')

library(data.table)
library(tidyverse)
library(ggrepel)

file_TF <- 'TF_names_v_1.01.txt'

TF <- fread(file_TF, header = FALSE, col.names = "TF_gene")

tmp <- sapply(TF$TF_gene, function(x) grep(x, PCO_trans_uniq$neargene) )
tmp_df <- tibble(TF_gene = names(tmp),
                 l = sapply(tmp, length) ) %>%
  arrange(desc(l))
View(tmp_df)

gwasResults <- readRDS('manhattan.rds')
PCO_trans_uniq <- fread("PCO.trans.uniq.txt")
don <- readRDS('don.rds')

anno_dat = PCO_trans_uniq %>% group_by(nearestgene) %>% summarise(SNP_annotate_neargene = SNP[which.min(pvalue)] )


sapply(c("WDR4", "S100P", "AC093323.1", "xxx"), function(x) grep(x, PCO_trans_uniq$neargene[c(1, 50)]) )

gene_of_interest <- c("NFKBIA", "PLAGL1", "NFE2", "IKZF1", "KLF14", "NFKB1", "ZNF229", "BAZ2B", #TF
                      "ARHGEF3", "SENP7" )


cis_tf <- apply(PCO_trans_uniq, 1, function(x){
  tmp_gene = strsplit(x["neargene"], split = ";")[[1]]
  ind = gene_of_interest %in% tmp_gene
  
  gene_of_interest[ind]
})
df_TF_label <- tibble("if_label" = as.logical(sapply(cis_tf, length )),
                      "TF_label" = NA )
df_TF_label$TF_label[df_TF_label$if_label] <- unlist(cis_tf)

PCO_trans_uniq <- cbind(PCO_trans_uniq, df_TF_label)


axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(chromosome) + min(chromosome) ) / 2 )
label = PCO_trans_uniq[match(don[don$is_annotate_neargene=="yes", "SNP"], PCO_trans_uniq$SNP), ]$nearestgene
label[label=="NA"] = ""

anno_dat = PCO_trans_uniq %>% group_by(TF_label) %>% summarise(SNP_annotate_TF_label = SNP[which.min(pvalue)] ) %>% ungroup()
anno_dat <- anno_dat %>% filter( !is.na(TF_label) )

don <- don %>% mutate(is_annotate_TF = ifelse(SNP %in% anno_dat$SNP_annotate_TF_label, "yes", "no"))
TF_label = PCO_trans_uniq[match(don[don$is_annotate_TF=="yes", "SNP"], PCO_trans_uniq$SNP), ]$TF_label
TF_label[is.na(TF_label)] = ""


manhattan_plot <- ggplot(don, aes(x=chromosome, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#bfbfbf", "#8c8c8c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(6.5, 21) ) +
  geom_point(data=subset(don, is_annotate_TF=="yes"),
             color = "#b20000", fill = "#b20000", shape=23, size = 3)

fig <- manhattan_plot +
  labs(x = "Chromosome", y = "-Log10(P)") +
  geom_text_repel( data=subset(don, is_annotate_TF=="yes"),
                   aes(label=TF_label),
                   segment.colour="black",
                   size = 4,
                   min.segment.length = 0,
                   max.overlaps = 1,
                   fill = NA,
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
  ) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour="black", size = 0.7),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text=element_text(colour = "black", size=14),
        axis.title.y = element_text(angle=90,vjust =2, size=16),
        axis.title.x = element_text(vjust = -0.2, size=16) )

saveRDS(fig, "fig3_manhanttan_plot.rds")
saveRDS(TF_label, "fig3_obj1.rds")
ggsave('fig3_manhanttan_plot.png', fig, width = 7, height = 4)

