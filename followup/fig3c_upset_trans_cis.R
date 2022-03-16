rm(list = ls())

library(data.table)
library(UpSetR)

file_trans <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/LD.prun.in.chr.module.perm10.fdr10.txt'
file_cis_e <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
file_cis_s <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_sQTL_signif_variant_gene_pairs.txt.gz'

trans <- fread(file_trans, header = FALSE)
cis_e <- fread(file_cis_e, header = TRUE)
cis_s <- fread(file_cis_s, header = TRUE)

listInput = list("trans" = trans$V1,
                 "cis_eQTLs" = unique(cis_e$sid),
                 "cis_sQTLs" = unique(cis_s$sid))
fig <- upset(fromList(listInput),
             order.by = "freq",
             intersections = list(list("trans"),
                                  list("trans", "cis_sQTLs"),
                                  list("trans", "cis_eQTLs"),
                                  list("trans", "cis_sQTLs", "cis_eQTLs")
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

saveRDS(fig, "upset.rds")

pdf(file="upset.pdf", width = 6, height = 5) # or other device
fig
dev.off()
