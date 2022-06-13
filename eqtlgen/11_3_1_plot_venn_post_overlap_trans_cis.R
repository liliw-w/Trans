###############################################################
########### Plot the overlapps of trans with cis ###########
########### by Venn diagram ###########
###############################################################
rm(list = ls())
library(data.table)
library(VennDiagram)


### sets
file_trans <- '/scratch/midway2/liliw1/eQTGen_est_Sigma/postanalysis/LD.prun.in.chr.module.txt'
file_cis_e <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
file_cis_s <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_sQTL_signif_variant_gene_pairs.txt.gz'


### read sets
trans <- fread(file_trans, header = FALSE)
cis_e <- fread(file_cis_e, header = TRUE)
cis_s <- fread(file_cis_s, header = TRUE)


### venn input
listInput = list("Trans" = trans$V1,
                 "Cis-eQTLs" = unique(cis_e$sid),
                 "Cis-sQTLs" = unique(cis_s$sid))


# Prepare a palette of 3 colors with R colorbrewer:
myCol <- RColorBrewer::brewer.pal(length(listInput), "Pastel2")


# Chart
venn_obj <- venn.diagram(
  x = listInput,
  category.names = c("Trans" , "Cis-eQTLs" , "Cis-sQTLs"),
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
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


### save figure
saveRDS(venn_obj, "plot/trans_cis_venn.rds")

pdf(file = "plot/trans_cis_venn.pdf", width = 2, height = 2)
grid::grid.draw(venn_obj)
dev.off()

