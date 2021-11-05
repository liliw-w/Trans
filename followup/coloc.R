rm(list = ls())
library(data.table)
library(tidyverse)
library(coloc)

########## files and parameter, read files ##########
file_snp_meta_qtl = "/scratch/midway2/liliw1/coloc/snp.meta.qtl.txt.gz"
file_qtlColocReg = "/scratch/midway2/liliw1/coloc/qtlColocReg.txt.gz"

dir_p = "/project2/xuanyao/llw/DGN_no_filter_on_mappability/p/"
Nmodule = 166
module_seq = 1:Nmodule
p_included_thre = 1e-5
regionDis = 1e5


########## Define regions ##########
qtlColocReg = NULL
for(module in module_seq){
  ### pvalue files names
  file.p = as.character(outer(module, 1:22, FUN = function(x, y) paste0("p.module", x, ".chr", y, ".rds")))

  ### read p
  p.obs = rbindlist(lapply(file.p, function(x)
  {tmp_y=readRDS(paste0(dir_p, x));
  print(x);
  if(!is.null(tmp_y)){
    as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)
  }
  }))
  # p.obs$V2[p.obs$V2==0] = 1e-20
  p.obs = p.obs %>% rename(Signal = V1, Pval = V2)

  if(min(p.obs$Pval) > p_included_thre) next

  ### reformat
  tmpqtlColocReg = p.obs %>%
    separate(col = Signal, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE) %>%
    unite(col = "SNP_ID", c("Chr", "Pos"), sep = ":", remove = FALSE) %>%
    arrange(Pval) %>%
    mutate(Region = NA) %>% mutate(Included = FALSE)
  # tmpqtlColocReg = tmpqtlColocReg %>% filter(Pval < p_included_thre)

  ### Find the lead SNP of a region and define the region
  lSignal = tmpqtlColocReg[1, ]
  while (lSignal$Pval <= p_included_thre) {
    indReg = tmpqtlColocReg$Chr == lSignal$Chr & abs(tmpqtlColocReg$Pos-lSignal$Pos) < regionDis/2
    tmpqtlColocReg[indReg, "Included"] = TRUE
    tmpqtlColocReg[indReg, "Region"] = lSignal$Signal

    qtlColocReg = rbind(qtlColocReg, tmpqtlColocReg[tmpqtlColocReg$Included, ])

    tmpqtlColocReg = tmpqtlColocReg[!tmpqtlColocReg$Included, ]
    lSignal = tmpqtlColocReg[1, ]

  }

  fwrite(qtlColocReg, file_qtlColocReg, quote = FALSE, sep = "\t")
  cat("coloc region defined for module", module, ", out of", Nmodule, "modules.", "\n")

}


## overlap with gwas, add rsid from gwas
snp_meta_qtl = fread(file_snp_meta_qtl)
qtlColocReg = left_join(x = qtlColocReg, y = snp_meta_qtl, by = c("SNP_ID", "Chr", "Pos"))

fwrite(qtlColocReg, file_qtlColocReg, quote = FALSE, sep = "\t")


## add MAF column, plink?


########## Summarize the region info ##########
### 0. # overlapped SNPs in both QTL and GWAS

### 1. #regions for each module?

### 2. #SNPs in each region



########## prepare coloc files and run coloc ##########
D1 = list("pvalues",
          "N",
          "MAF",
          "type",
          "snp")
D2 = list("pvalues",
          "N",
          "MAF",
          "type",
          "snp")
maf

check_dataset(D1)
check_dataset(D2)

## coloc
coloc_res = coloc.abf(D1, D2, MAF = maf)

## follow up analysis

## save results

