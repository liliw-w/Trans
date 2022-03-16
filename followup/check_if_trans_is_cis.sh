######## BASH ###########
######## Check if trans-eQTLs are cis-eQTLs ###########
file_trans=/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/LD.prun.in.chr.module.perm10.fdr10.txt
file_cis=/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_sQTL_signif_variant_gene_pairs.txt.gz
prefix='sQTL'
if_GTEx_cis='no'

# cd /scratch/midway2/liliw1/trans_cis/

if [[ -z  ${if_GTEx_cis} ]] || [[ -z  ${file_trans} ]] || [[ -z  ${file_cis} ]] || [[ -z ${prefix} ]]
then
  prefix='sQTL'
  if_GTEx_cis='no'
  file_cis=/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_sQTL_signif_variant_gene_pairs.txt.gz
  #/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz
  #/project2/xuanyao/data/GTEx_v7/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz
  file_trans=/project2/xuanyao/llw/DGN_PCO.lambda.01_real/postanalysis/LD.prun.in.chr.module.perm10.fdr10.txt
  #/project2/xuanyao/llw/eQTLGen_DGN_PCO.lambda.01/postanalysis/LD.prun.in.chr.module.perm10.txt
  #/project2/xuanyao/llw/DGN_PCO.lambda.01/postanalysis/LD.prun.in.chr.module.perm10.txt
fi

if [[ ${if_GTEx_cis} = 'no' ]]
then
  zcat ${file_cis} | grep -f ${file_trans} | sed "1 i $(zcat ${file_cis} | head -n1)" > trans_cis_${prefix}.txt
  cat trans_cis_${prefix}.txt | tail -n+2 | cut -f 2 | sort | uniq > signals_${prefix}.txt
  cat trans_cis_${prefix}.txt | tail -n+2 | cut -f 1 | sort | uniq > Genes_${prefix}.txt

  Ntrans=$(cat ${file_trans} | wc -l)
  Ntrans_cis=$(cat signals_${prefix}.txt | wc -l)
  Ncis_eGenes=$(cat Genes_${prefix}.txt | wc -l)
  cis_eGenes=$(cat Genes_${prefix}.txt)
else
  file_anno=/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt
  #/project2/xuanyao/data/GTEx_v7/Whole_Blood.v7.egenes.txt.gz

  cat ${file_trans} | sed 's\:\_\g' > tmp.txt
  zcat ${file_cis} | grep -f tmp.txt | sed "1 i $(zcat ${file_cis} | head -n1)" > trans_cis_${prefix}.txt
  cat trans_cis_${prefix}.txt | tail -n+2 | cut -f 1 | sort | uniq > signals_${prefix}.txt
  cat trans_cis_${prefix}.txt | tail -n+2 | cut -f 2 | sort | uniq > Genes_${prefix}.txt
  cat ${file_anno} | grep -w -f Genes_${prefix}.txt | sed "1 i $(head -n1 ${file_anno})" > Genes_${prefix}_annoted.txt
  rm tmp.txt

  Ntrans=$(cat ${file_trans} | wc -l)
  Ntrans_cis=$(cat signals_${prefix}.txt | wc -l)
  Ncis_eGenes=$(cat Genes_${prefix}.txt | wc -l)
  cis_eGenes=$(cat Genes_${prefix}_annoted.txt | tail -n+2 | cut -f 2)
fi
echo "${Ntrans_cis} trans-eQTLs (out of ${Ntrans}) are also cis-${prefix}, corresponding to ${Ncis_eGenes} cis-${prefix}-Genes. \n"; \
echo "cis-${prefix}-Genes include: \n ${cis_eGenes}"
