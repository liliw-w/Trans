##############################################
########### define trans-eQTLs regions for coloc ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
script_cal_maf <- "cal_maf.sh"
dir_geno <- "/project2/xuanyao/llw/DGN/data/"
dir_p <- "/project2/xuanyao/llw/MODULES/MSigDB/p/"
Nmodule <- 50
module_seq <- 1:Nmodule
p_included_thre <- 1e-5
dis_reg <- 1e5

file_snp_meta_qtl <- "/project2/xuanyao/llw/coloc/snp.meta.qtl.txt.gz"
file_qtlMaf <- "/project2/xuanyao/llw/DGN/data/chr_all.frq"

## output -----
file_qtlColocReg <- "/scratch/midway2/liliw1/coloc_MSigDB/qtlColocReg.txt.gz"


# Define regions -----
qtlColocReg <- NULL
for(module in module_seq){
  ## read p -----
  file.p <- as.character(outer(module, 1:22, FUN = function(x, y) paste0("p.module", x, ".chr", y, ".rds")))
  
  p.obs <- rbindlist(lapply(file.p, function(x)
  {tmp_y=readRDS(paste0(dir_p, x));
  print(x);
  if(!is.null(tmp_y)){
    as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)
  }
  }))
  # p.obs$V2[p.obs$V2==0] <- 1e-20
  p.obs <- p.obs %>% rename(Signal = V1, Pval = V2)

  if(min(p.obs$Pval) > p_included_thre) next

  ## reformat signals and snps -----
  tmpqtlColocReg <- p.obs %>%
    separate(col = Signal, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE) %>%
    unite(col = "SNP_ID", c("Chr", "Pos"), sep = ":", remove = FALSE) %>%
    arrange(Pval) %>%
    mutate(Region = "") %>% mutate(Included = FALSE)
  # tmpqtlColocReg <- tmpqtlColocReg %>% filter(Pval < p_included_thre)

  ## Find the lead SNP of a region and define the region -----
  lSignal <- tmpqtlColocReg[1, ]
  while (lSignal$Pval <= p_included_thre) {
    indReg <- tmpqtlColocReg$Chr == lSignal$Chr & abs(tmpqtlColocReg$Pos-lSignal$Pos) < dis_reg/2
    tmpqtlColocReg[indReg, "Included"] <- TRUE
    tmpqtlColocReg[indReg, "Region"] <- lSignal$Signal

    qtlColocReg <- rbind(qtlColocReg, tmpqtlColocReg[tmpqtlColocReg$Included, ])

    tmpqtlColocReg <- tmpqtlColocReg[!tmpqtlColocReg$Included, ]
    lSignal <- tmpqtlColocReg[1, ]

  }
  
  #fwrite(qtlColocReg, file_qtlColocReg, quote = FALSE, sep = "\t")
  cat("coloc region defined for module", module, ", out of", Nmodule, "modules.", "\n")

}


# overlap with gwas snps, add rsid from gwas -----
snp_meta_qtl <- fread(file_snp_meta_qtl)
qtlColocReg <- left_join(x = qtlColocReg, y = snp_meta_qtl, by = c("SNP_ID", "Chr", "Pos"))


# add MAF column, MAF calcualted from plink -----
if(!file.exists(file_qtlMaf)){
  cmd <- paste("bash", script_cal_maf, dir_geno)
  system(cmd)
}

qtlMaf <- fread(file_qtlMaf, header = TRUE)
qtlMaf <- qtlMaf[!duplicated(qtlMaf$SNP), ]
if( !all(qtlColocReg$SNP_ID %in% qtlMaf$SNP) ) stop("Not all SNPs have maf!")

qtlColocReg <- qtlColocReg %>%
  left_join(y = qtlMaf, by = c("SNP_ID" = "SNP", "Chr" = "CHR")) %>%
  mutate(A1 = NULL, A2 = NULL, NCHROBS = NULL) %>%
  arrange(Pval)


# print out key message or write out -----
fwrite(qtlColocReg, file_qtlColocReg, quote = FALSE, sep = "\t")

