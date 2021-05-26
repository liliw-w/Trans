######## Find all cis-eQTLs for every eGene ###########
require(data.table)

args = commandArgs(trailingOnly = TRUE)
file_input = args[1]
file_all = args[2]
file_eGene = args[3]
file_eQTL = args[4]
fdr = args[5]

if(is.na(file_input) | is.na(file_all) | is.na(file_eGene) | is.na(file_eQTL) | is.na(fdr)){
  file_input = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_eQTL.txt.gz'
  file_all = '~/xuanyao_llw/DGN_PCO.lambda.01/nom_20phenoPC_allChr_eQTL.txt.gz'
  file_eGene = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_egenes.txt.gz'
  file_eQTL = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
  fdr = 0.05
}


## Append pval_nominal_threshold to each gene
eGene_orig = fread(file_input, header=TRUE)
eGene = eGene_orig[eGene_orig$qval < fdr, ]

pt = mean(c(max(eGene$bpval), min(eGene$bpval)))
eGene$pval_nominal_threshold = signif(qbeta(pt,
                                            eGene$shape1, eGene$shape2, ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
fwrite(eGene, file_eGene, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


## Find cis-eQTLs for each gene
assoc_all = fread(file_all, header = FALSE, sep=" ")
setnames(assoc_all, c("pid", "sid", "dist", "npval", "slope"))
print(head(assoc_all))

sig_egene = unique(eGene$pid)
assoc_all = assoc_all[assoc_all$pid %in% sig_egene, ]
assoc_all$pval_nominal_threshold = eGene[match(assoc_all$pid, eGene$pid), pval_nominal_threshold]
assoc_all$if_QTL = assoc_all$npval <= assoc_all$pval_nominal_threshold

fwrite(assoc_all[assoc_all$if_QTL, ], file_eQTL, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
