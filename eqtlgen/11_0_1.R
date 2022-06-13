###############################################################
########### Plot the overlapps of trans with DGN trans ###########
########### by Venn diagram ###########
###############################################################
rm(list = ls())
library(VennDiagram)
library(grid)


### sets
file_trans <- '/scratch/midway2/liliw1/eQTGen_est_Sigma/postanalysis/LD.prun.in.chr.module.txt'
file_trans2 <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/LD.prun.in.chr.module.perm10.fdr10.txt'


### read sets
trans <- fread(file_trans, header = FALSE)
trans2 <- fread(file_trans2, header = FALSE)


### venn input
listInput = list("Trans-eQTLGen" = unique(trans$V1),
                 "Trans-DGN" = unique(trans2$V1))


# Prepare a palette of 3 colors with R colorbrewer:
myCol <- RColorBrewer::brewer.pal(8, "Pastel2")


# Chart
venn_obj <- venn.diagram(
  x = listInput,
  category.names = c("Trans-eQTLGen", "Trans-DGN"),
  filename = NULL,
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 200, 
  width = 200, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:length(listInput)],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #rotation = 1
)


### save figure
saveRDS(venn_obj, "plot/trans_trans_venn.rds")

pdf(file = "plot/trans_trans_venn.pdf", width = 2, height = 2)
grid.draw(venn_obj)
dev.off()

