rm(list = ls())
library(data.table)
library(tidyverse)


dir_cis_qtl <- "/project2/xuanyao/llw/DGN_sQTL"
file_intron_to_gene <- "/project2/xuanyao/llw/DGN_sQTL/merged_all_intron_to_gene.txt.gz"
#file_cis_qtl <- "/project2/xuanyao/llw/DGN_PCO.lambda.01/nom_10phenoPC_allChr.sQTL.txt.gz"



### trans- regions
#### sorted by chr



### cis- qtl regions
file_cis_qtl_list <- list.files(path = dir_cis_qtl,
                            pattern = "^nom_10phenoPC_.*txt.gz$",
                            full.names = TRUE)
intron_to_gene <- fread(file_intron_to_gene,
                        col.names = c("chr_gene", "start_gene", "end_gene", "gene",
                                      "chr_intron", "start_intron", "end_intron", "intron",
                                      "nSNP") )

#### read by chr
for(chr in 1:22){
  file_cis_qtl_chr <- grep(paste0("nom_10phenoPC_", chr, ".maf5e-2.txt.gz"),
                           file_cis_qtl_list,
                           value = TRUE)
  
  cis_qtl_chr <- fread(file_cis_qtl_chr,
                       col.names = c("pid", "sid", "dist", "npval", "slope"))
  
  ### Add gene name to intron
  intron_list <- distinct(cis_qtl_chr, pid) %>%
    separate("pid",
             c("chr", "start", "end", NA, NA),
             sep = ":", remove = FALSE, convert = TRUE) %>%
    unite("intron", chr, start, end, sep = ":")
  intron_list <- intron_list %>%
    left_join(intron_to_gene, by = "intron")
  
  fwrite(intron_list, file = paste0("intron_to_gene_chr", chr, ".txt.gz") )
  
  cat("CHR", chr, "\n")
}


### two-way overlap


### remove reg with less than ??? SNPs


### light blue regions (cis- regions)


### light green regions, candidate coloc regions



