##############################################
########## Prepare coloc regions for gwas traits ##########
########## These traits are part of the 72 traits ##########
########## Focus on the 14 autoimmune diseases ##########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# Trait: RA, add beta, se columns -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid24390342_RA_GWASmeta_European_v2.txt.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## beta
gwas$beta <- log(gwas[["OR(A1)"]])

## se
gwas$se <- with(gwas, ( log(gwas[["OR_95%CIup"]]) - log(gwas[["OR_95%CIlow"]]) )/3.92 )

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")




# Trait: IIBDGC - IBD, add beta -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid26192919_EUR.IBD.gwas_info03_filtered.assoc.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## beta
gwas$beta <- log(gwas[["OR"]])

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")


# Trait: IIBDGC - CD, add beta -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid26192919_EUR.CD.gwas_info03_filtered.assoc.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## beta
gwas$beta <- log(gwas[["OR"]])

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")


# Trait: IIBDGC - UC, add beta -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid26192919_EUR.UC.gwas_info03_filtered.assoc.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## beta
gwas$beta <- log(gwas[["OR"]])

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")



# Trait: Lange IBD, add chr and pos -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid28067908_ibd_build37_59957_20161107.txt.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## chr and pos
gwas <- gwas %>%
  separate(col = MarkerName, into = "SNP", sep = "_", remove = FALSE, convert = TRUE) %>%
  separate(col = SNP, into = c("Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")



# Trait: Lange CD, add chr and pos -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid28067908_cd_build37_40266_20161107.txt.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## chr and pos
gwas <- gwas %>%
  separate(col = MarkerName, into = "SNP", sep = "_", remove = FALSE, convert = TRUE) %>%
  separate(col = SNP, into = c("Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")



# Trait: Lange UC, add chr and pos -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid28067908_uc_build37_45975_20161107.txt.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## chr and pos
gwas <- gwas %>%
  separate(col = MarkerName, into = "SNP", sep = "_", remove = FALSE, convert = TRUE) %>%
  separate(col = SNP, into = c("Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")



# Trait: MS, add beta, se columns -----
file_gwas <- "/project2/xuanyao/llw/GWAS/pmid31604244_MS_15_discovery_metav3.0.meta.gz"

## original gwas file
gwas <- fread(file = file_gwas)

## beta
gwas$beta <- log(gwas[["OR"]])

## se
gwas$abs_z <- with(gwas, sqrt(qchisq(1 - P, 1)))
gwas$se <- with(gwas, abs(beta / abs_z) )
gwas$abs_z <- NULL

## save new formatted gwas file
fwrite(gwas, file_gwas, quote = FALSE, sep = "\t")

