###############################################################
########### Plot the overlapps of pc1 signals and trans-PCO signals ###########
########### by upset plot ###########
###############################################################
rm(list = ls())

library(data.table)
library(UpSetR)


### upset input sets -----
file_pco_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_pc1_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability_PC1/FDR/signals.chr.module.perm10.txt'

file_plot_upset_dat <- 'pc1/upset.rds'
file_plot_upset <- 'pc1/upset.pdf'


# read files -----
pco_sig <- fread(file_pco_sig, header = FALSE, col.names = c('signal', 'p', 'q'))
pc1_sig <- fread(file_pc1_sig, header = FALSE, col.names = c('signal', 'p', 'q'))


# data sets for upset plot -----
listInput <- list(
  "Trans_PCO" = pull(pco_sig, signal),
  "PC1" = pull(pc1_sig, signal)
)


# plot -----
fig <- upset(fromList(listInput),
             order.by = "freq",
             intersections = list(list("Trans_PCO", "PC1"),
                                  list("Trans_PCO"),
                                  list("PC1")
             ),
             decreasing = TRUE,
             mainbar.y.label = "Shared Signal",
             sets.x.label = "Signal",
             point.size = 12,
             line.size = 2.5,
             mb.ratio = c(0.6, 0.4),
             text.scale = c(2, 2, 2, 1, 2, 2),
             set_size.show = TRUE,
             #set_size.angles = 30,
             set_size.numbers_size = 5)

# save plot -----
saveRDS(fig, file_plot_upset_dat)

pdf(file = file_plot_upset, height = 4, width = 5)
fig
dev.off()

