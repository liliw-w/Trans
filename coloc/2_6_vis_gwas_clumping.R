##############################################
########### merge ld clumped files and visualize ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)

source('/project2/xuanyao/llw/Trans/plot/theme_my_pub.R')
theme_set(
  theme_my_pub(legend.position = "right")
)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

gwas_trait_type <- "continuous"

file_pheno_manifest <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv"

## output -----


# read files -----
pheno_manifest <- data.table::fread(file_pheno_manifest)
gwasPhenocode <- pheno_manifest$`GWAS ID`
file_clump_res <- str_glue('ld_clump/clumped_snp_allchr_{gwas_trait_type}-{gwasPhenocode}.txt')


# vis -----
## data for vis -----
plt_dat <- lapply(
  file_clump_res,
  function(x) {
    data.table::fread(x, header = FALSE) %>%
      nrow() %>%
      enframe(name = NULL, value = "n_ld_clumped")
  }
) %>%
  setNames(gwasPhenocode) %>%
  bind_rows(.id = "GWAS ID") %>%
  mutate(`GWAS ID` = as.numeric(`GWAS ID`)) %>%
  left_join(
    pheno_manifest,
    by = "GWAS ID"
  )

## vis -----
ggplot(plt_dat, aes(y = `Trait Abbreviation`, x = n_ld_clumped)) +
  geom_col() +
  geom_text(aes(label = n_ld_clumped), nudge_x = 200, size = 2) +
  labs(x = "Number of LD clumped GWAS loci", y = NULL)


## save fig -----
ggsave(
  "ld_clump/num_ld_clumped_gwas_loci.pdf",
  width = 2, height = 3
)

