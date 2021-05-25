######## BASH ###########
######## Check if trans-eQTLs are cis-eQTLs ###########
file_trans=$1
file_cis=$2
if_GTEx_cis=$3

if [[ -z  ${if_GTEx_cis} ]] || [[ -z  ${file_trans} ]] || [[ -z  ${file_cis} ]]
then
  if_GTEx_cis='no'
  file_trans=/project2/xuanyao/llw/DGN_PCO.lambda.01/postanalysis/LD.prun.in.chr.module.perm10.txt
  #/project2/xuanyao/llw/eQTLGen_DGN_PCO.lambda.01/postanalysis/LD.prun.in.chr.module.perm10.txt
  #/project2/xuanyao/llw/DGN_PCO.lambda.01/postanalysis/LD.prun.in.chr.module.perm10.txt
  file_cis=/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz
  #/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz
  #/project2/xuanyao/data/GTEx_v7/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz
fi

if [[ ${if_GTEx_cis} = 'no' ]]
then
  zcat ${file_cis} | grep -f ${file_trans} | sed "1 i $(zcat ${file_cis} | head -n1)" > trans_cis.txt
  cat trans_cis.txt | tail -n+2 | cut -f 2 | sort | uniq > cis_transeQTLs.txt
  cat trans_cis.txt | tail -n+2 | cut -f 1 | sort | uniq > cis_eGenes.txt

  Ntrans=$(cat ${file_trans} | wc -l)
  Ntrans_cis=$(cat cis_transeQTLs.txt | wc -l)
  Ncis_eGenes=$(cat cis_eGenes.txt | wc -l)
  cis_eGenes=$(cat cis_eGenes.txt)
else
  file_anno=/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt
  #/project2/xuanyao/data/GTEx_v7/Whole_Blood.v7.egenes.txt.gz

  cat ${file_trans} | sed 's\:\_\g' > tmp.txt
  zcat ${file_cis} | grep -f tmp.txt | sed "1 i $(zcat ${file_cis} | head -n1)" > trans_cis.txt
  cat trans_cis.txt | tail -n+2 | cut -f 1 | sort | uniq > cis_transeQTLs.txt
  cat trans_cis.txt | tail -n+2 | cut -f 2 | sort | uniq > cis_eGenes.txt
  cat ${file_anno} | grep -w -f cis_eGenes.txt | sed "1 i $(head -n1 ${file_anno})" > cis_eGenes_annoted.txt
  rm tmp.txt

  Ntrans=$(cat ${file_trans} | wc -l)
  Ntrans_cis=$(cat cis_transeQTLs.txt | wc -l)
  Ncis_eGenes=$(cat cis_eGenes.txt | wc -l)
  cis_eGenes=$(cat cis_eGenes_annoted.txt | tail -n+2 | cut -f 2)
fi
echo "${Ntrans_cis} trans-eQTLs (out of ${Ntrans}) are also cis-eQTLs, corresponding to ${Ncis_eGenes} cis-eGenes. \n"; \
echo "cis-eGenes include: \n ${cis_eGenes}"
