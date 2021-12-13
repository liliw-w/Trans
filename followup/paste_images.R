library(png)

gwasPhenocode = "30090"
tmp_file = list.files("~/Downloads", "^pheno30090.*.reg.coloc.png$")

for(tmp in tmp_file){
  tmp_split = strsplit(tmp, ".", fixed = TRUE)[[1]]
  reg = tmp_split[2]

  file_coloc = paste0("~/Downloads/pheno", gwasPhenocode, ".", reg, ".reg.coloc.png")
  file_mr = paste0("pheno", gwasPhenocode, ".", reg, ".mr.png")

  img_coloc <- readPNG(file_coloc)
  img_mr <- readPNG(file_mr)

  pdf(paste0("pheno", gwasPhenocode, ".", reg, ".pdf"))
  plot(NA, xlim = c(0, 8), ylim = c(0, 2), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  rasterImage(img_coloc, 0, 0, 4, 1.5)
  rasterImage(img_mr, 4, 0, 8, 1.5)
  dev.off()

  cat("Region", reg, "is done. \n")
}



