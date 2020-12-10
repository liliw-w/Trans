#!/bin/bash 

dir_cov=/project2/xuanyao/data/GTEx_V8/expression/GTEx_Analysis_v8_eQTL_covariates
dir_exp=/project2/xuanyao/data/GTEx_V8/expression/GTEx_Analysis_v8_eQTL_expression_matrices
dir_GTEx=/project2/xuanyao/llw/GTEx_v8

all_tissue=$(ls ${dir_cov} | cut -d "." -f 1)
for tissue in ${all_tissue}
do
dir_data=${dir_GTEx}/${tissue}
prefix=${tissue}.v8
echo ${dir_data}

if ! [ -d "${dir_data}/data" ]; then
    if ! [ -d "${dir_data}" ]; then mkdir ${dir_data}; fi
    mkdir ${dir_data}/data
fi

cd ${dir_data}/data
cp ${dir_cov}/${prefix}* ./covariates.txt
cp ${dir_exp}/${prefix}* .

# write gene meta file
zcat ${prefix}*.bed.gz | \
cut -f 1-4 | tail -n+2 | \
awk 'BEGIN {printf("#chr\tstart\tend\tgene\n")} {print $0}' > \
gene.meta.txt

Rscript ${dir_GTEx}/prep.exp.R ${prefix}*.bed.gz 'ex.rdata'

rm ${prefix}*

cp /project2/xuanyao/llw/TCGA/data/pseudogenes.txt .
cp /project2/xuanyao/llw/TCGA/data/mappability.txt .
cp /project2/xuanyao/llw/TCGA/data/cross.mappable.genes.rds .

done
