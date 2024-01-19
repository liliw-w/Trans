# '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'
# '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/qtl_table.txt'


file_qtl <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'

file_unique_snp <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/ld_clump/uniq_snp.txt'
file_with_p <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/ld_clump/ld_clump_esnp.txt'


qtl <- data.table::fread(file_qtl)

# unique snps
distinct(qtl, SNP) %>% 
  data.table::fwrite(
    ., 
    file = file_unique_snp,
    sep = '\t', quote = FALSE, col.name = FALSE
  )


# snp/snp pair with p
group_by(qtl, SNP) %>% 
  summarise(
    P = min(P)
  ) %>% 
  arrange(SNP) %>% 
  data.table::fwrite(
    ., 
    file = file_with_p,
    sep = '\t', quote = FALSE
  )



