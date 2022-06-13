##############################################
########### Look at the replication of ###########
########### eQTLGen trans signal pairs in ###########
########### DGN trans signals ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
fdr_level <- 0.1
ratio <- 50

file_eqtlgen_sig <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')

file_dgn_all_snp <- "/project2/xuanyao/llw/eQTLGen_DGN/DGN.all_snp.txt"
file_p_all <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/q.chr.module.perm2.rds'

file_rep <- 'postanalysis/rep_eqtlgen_in_dgn.txt'


# read files -----
eqtlgen_sig <- fread(file_eqtlgen_sig, header = TRUE)
dgn_all_snp <- fread(file_dgn_all_snp, header = FALSE)
p_all <- readRDS(file_p_all)



# organize data -----
eqtlgen_sig$module <- paste0("module", eqtlgen_sig$module)

dgn_all_snp_uniq <- dgn_all_snp$V1

p_all <- p_all %>%
  separate(snp,
           into = c('module', 'chr', 'pos'),
           sep = ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_id", c(chr, pos), sep = ":", remove = FALSE)


# filter eqtlgen signals that were also anlyzed by dgn  -----
# and add dgn p values
eqtlgen_sig_dgn <- eqtlgen_sig %>%
  filter(meta %in% dgn_all_snp_uniq) %>%
  left_join(p_all %>% select(module, p, SNP_id),
            by = c("module", "meta" = "SNP_id"),
            suffix = c(".eqtlgen", ".dgn")) %>%
  replace_na(list(p.dgn = 1)) %>%
  arrange(desc(p.dgn))


# define new p threshold for replicated signals  -----
p_thre <- fdr_level/nrow(eqtlgen_sig_dgn)


eqtlgen_sig_dgn <- eqtlgen_sig_dgn %>%
  mutate("if_rep" = p.dgn < p_thre) %>%
  arrange(desc(if_rep), SNPChr, SNPPos, p.eqtlgen, p.dgn) %>%
  relocate(p.eqtlgen, .before = p.dgn)



# print out key message and write out  -----
cat(
  "Out of", nrow(eqtlgen_sig), "eQTLGen (SNP, module) signal pairs,", nrow(eqtlgen_sig_dgn), "of them are also analyzed in DGN. \n\n",
  
  "Based on these overlaps, the updated p-values threshold for defining replication is:", p_thre, ", under FDR level", fdr_level, ". \n\n",
  
  "As a result, out of", nrow(eqtlgen_sig_dgn), "eQTLGen signals,", sum(eqtlgen_sig_dgn$if_rep), "are replicated in DGN. \n\n",
  
  "The p-values of these signals in DGN range: (", min(eqtlgen_sig_dgn[eqtlgen_sig_dgn$if_rep, "p.dgn"]), ",", max(eqtlgen_sig_dgn[eqtlgen_sig_dgn$if_rep, "p.dgn"]), "). \n\n"
)


fwrite(eqtlgen_sig_dgn, file_rep,
       sep = "\t", quote = FALSE)


