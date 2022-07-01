###############################################################
########### Plot the overlapps of trans with cis ###########
########### by upset plot ###########
###############################################################
rm(list = ls())

library(data.table)
library(UpSetR)


### sets
file_trans <- '/project2/xuanyao/llw/eQTLGen_est_Sigma/postanalysis/LD.prun.in.chr.module.txt'
file_cis_e <- '/project2/xuanyao/llw/DGN_sumstat/DGN_eQTL/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
file_cis_s <- '/project2/xuanyao/llw/DGN_sumstat/DGN_sQTL/DGN_sQTL_signif_variant_gene_pairs.txt.gz'

### read sets
trans <- fread(file_trans, header = FALSE)
cis_e <- fread(file_cis_e, header = TRUE)
cis_s <- fread(file_cis_s, header = TRUE)


### upset input
listInput <- list("Trans" = trans$V1,
                 "Cis_eQTLs" = unique(cis_e$sid),
                 "Cis_sQTLs" = unique(cis_s$sid))

### plot
fig <- upset(fromList(listInput),
             order.by = "freq",
             intersections = list(list("Trans"),
                                  list("Trans", "Cis_sQTLs"),
                                  list("Trans", "Cis_eQTLs"),
                                  list("Trans", "Cis_sQTLs", "Cis_eQTLs")
             ),
             decreasing = TRUE,
             mainbar.y.label = "Shared QTL",
             sets.x.label = "QTL tested",
             point.size = 12,
             line.size = 2.5,
             mb.ratio = c(0.6, 0.4),
             text.scale = c(2, 2, 2, 1, 2, 2),
             set_size.show = TRUE,
             #set_size.angles = 30,
             set_size.numbers_size = 5)
fig


### save plot
saveRDS(fig, "plot/trans_cis_upset.rds")

pdf(file = "plot/trans_cis_upset.pdf", width = 6, height = 5) # or other device
fig
dev.off()

