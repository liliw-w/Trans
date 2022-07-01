###############################################################
########### Plot the overlapps of trans with cis ###########
########### by upset plot ###########
###############################################################
rm(list = ls())

library(data.table)
library(UpSetR)


### upset input sets
file_upset_inp <- '/project2/xuanyao/llw/coloc/cis/coloc_upset_set.rds'


file_plot_upset_dat <- '/project2/xuanyao/llw/coloc/cis/coloc_trans_cis_upset.rds'
file_plot_upset <- '/project2/xuanyao/llw/coloc/cis/coloc_trans_cis_upset.pdf'


# read file
listInput <- readRDS(file_upset_inp)

names(listInput) <- c("Coloc", "Cis_eQTLs", "Cis_sQTLs", "Trans_loci")


### plot
fig <- upset(fromList(listInput),
             order.by = "freq",
             intersections = list(list("Trans_loci"),
                                  list("Cis_eQTLs", "Trans_loci"),
                                  list("Cis_sQTLs", "Trans_loci"),
                                  list("Cis_eQTLs", "Cis_sQTLs", "Trans_loci")
             ),
             decreasing = TRUE,
             mainbar.y.label = "Shared Region",
             sets.x.label = "Region",
             point.size = 12,
             line.size = 2.5,
             mb.ratio = c(0.6, 0.4),
             text.scale = c(2, 2, 2, 1, 2, 2),
             set_size.show = TRUE,
             #set_size.angles = 30,
             set_size.numbers_size = 5)
fig


### save plot
saveRDS(fig, file_plot_upset_dat)

pdf(file = file_plot_upset, width = 6, height = 5) # or other device
fig
dev.off()

