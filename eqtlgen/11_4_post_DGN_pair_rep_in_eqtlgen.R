##############################################
########### Look at the replication of ###########
########### DGN trans signal pairs in ###########
########### eQTLGen trans signals ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
fdr_level <- 0.1
ratio <- 50

file_dgn_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

file_module_use <- paste0('postanalysis/module_use_ratio_', ratio, '.txt')
file_p_all <- 'p/p.module_all.Sigma_nullz.rds'

file_rep <- 'postanalysis/rep_dgn_in_eqtlgen.txt'


# read files -----
dgn_sig <- fread(file_dgn_sig, header = FALSE, col.names = c("signal", "p", "q") )

p_all <- readRDS(file_p_all)
eqtlgen_all_snp <- p_all %>% distinct(meta) %>% pull()
module_use <- fread(file_module_use, header = TRUE)


# organize data -----
p_all$module <- paste0("module", p_all$module)

dgn_sig <- dgn_sig %>%
  separate(signal,
           into = c('module', 'chr', 'pos'),
           sep = ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_id", c(chr, pos), sep = ":", remove = FALSE)



# filter dgn signals that were also anlyzed by eqtlgen  -----
# and add eqtlgen p values
dgn_sig_eqtlgen <- dgn_sig %>%
  filter(module %in% paste0("module", module_use$module) &
           SNP_id %in% eqtlgen_all_snp ) %>%
  left_join(p_all %>% select(module, p, meta),
            by = c("module", "SNP_id" = "meta"),
            suffix = c(".dgn", ".eqtlgen")) %>%
  arrange(desc(p.eqtlgen))



# define new p threshold for replicated signals  -----
p_thre <- fdr_level/nrow(dgn_sig_eqtlgen)

dgn_sig_rep_in_eqtlgen <- dgn_sig_eqtlgen %>% filter(p.eqtlgen < p_thre)


dgn_sig_eqtlgen <- dgn_sig_eqtlgen %>%
  mutate("if_rep" = p.eqtlgen < p_thre) %>%
  arrange(desc(if_rep), p.dgn, p.eqtlgen) %>%
  relocate(p.dgn, .before = p.eqtlgen)


# print out key message and write out  -----
cat(
  "Out of", nrow(dgn_sig), "DGN (SNP, module) signal pairs,", nrow(dgn_sig_eqtlgen), "of them are also analyzed in eQTLGen. \n\n",
  
  "Based on these overlaps, the updated p-values threshold for defining replication is:", p_thre, ", under FDR level", fdr_level, ". \n\n",
  
  "As a result, out of", nrow(dgn_sig_eqtlgen), "DGN signals,", sum(dgn_sig_eqtlgen$if_rep), "are replicated in eQTLGen. \n\n",
  
  "The p-values of these signals in eQTLGen range: (", min(dgn_sig_eqtlgen[dgn_sig_eqtlgen$if_rep, "p.eqtlgen"]), ",", max(dgn_sig_eqtlgen[dgn_sig_eqtlgen$if_rep, "p.eqtlgen"]), "). \n\n"
)

fwrite(dgn_sig_eqtlgen, file_rep,
       sep = "\t", quote = FALSE)


