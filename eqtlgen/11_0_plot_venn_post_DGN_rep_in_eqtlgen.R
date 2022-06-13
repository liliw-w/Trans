###############################################################
########### Plot the overlapps of trans with cis ###########
########### by Venn diagram ###########
###############################################################
rm(list = ls())
library(VennDiagram)


### venn input
listInput = list(
  "DGN Signal" = dgn_sig_uniq,
  "eQTLGen Signal" = eqtlgen_old_signal_uniq,
  "eQTLGen Rep" = eqtlgen_sig_uniq,
  
  "eQTLGen All" = eqtlgen_all_snp_uniq,
  "DGN All" = dgn_all_snp_uniq
)



# Prepare a palette of 3 colors with R colorbrewer:
myCol <- RColorBrewer::brewer.pal(length(listInput), "Pastel2")
myCol <- c("#ffc9c9", "#cddcf0", "#adc6e6")


# Chart
listInput1 <- listInput[c(1, 3, 2)]
venn_obj <- venn.diagram(
  x = listInput1,
  #category.names = names(listInput),
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
  fill = myCol,
  
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.35,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0),
  cat.dist = c(0.05, 0.04, 0.02),
  cat.fontfamily = "sans",
  
  margin = 0
)

pdf(file = "venn_signal.pdf", width = 1.5, height = 1)
grid::grid.draw(venn_obj)
dev.off()




# Chart
myCol <- c("#ffc9c9", "#e5b4b4", "#cddcf0")

listInput2 <- listInput[c(5, 1, 4)]
venn_obj <- venn.diagram(
  x = listInput2,
  category.names = names(listInput2),
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
  fill = myCol,
  
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.35,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0),
  cat.dist = c(0.05, 0.04, 0.02),
  cat.fontfamily = "sans",
  
  rotation = 90,
  
  margin = 0
)

pdf(file = "venn_all_dgn.pdf", width = 1.5, height = 1)
grid::grid.draw(venn_obj)
dev.off()




# Chart
myCol <- c("#ffc9c9", "#cddcf0", "#d7e3f3")

listInput3 <- listInput[c(5, 2, 4)]
venn_obj <- venn.diagram(
  x = listInput3,
  category.names = names(listInput3),
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
  fill = myCol,
  
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.35,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0),
  cat.dist = c(0.05, 0.04, 0.02),
  cat.fontfamily = "sans",
  
  rotation = 90,
  
  margin = 0
)

pdf(file = "venn_all_etlgen.pdf", width = 1.5, height = 1)
grid::grid.draw(venn_obj)
dev.off()


