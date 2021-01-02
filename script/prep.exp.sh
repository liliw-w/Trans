#!/bin/bash 

dir_cov=/project2/xuanyao/data/GTEx_v8/expression/GTEx_Analysis_v8_eQTL_covariates
dir_exp=/project2/xuanyao/data/GTEx_v8/expression/GTEx_Analysis_v8_eQTL_expression_matrices
dir_GTEx=/scratch/midway2/liliw1/GTEx_v8/

all_tissue=$(ls ${dir_cov} | cut -d "." -f 1)
for tissue in ${all_tissue}
do
dir_data=${dir_GTEx}${tissue}
prefix=${tissue}.v8
echo ${dir_data}

if ! [ -d "${dir_data}/data" ]; then
    if ! [ -d "${dir_data}" ]; then mkdir ${dir_data}; fi
    mkdir ${dir_data}/data
fi

cp ${dir_cov}/${prefix}* ./data/covariates.txt
cp ${dir_exp}/${prefix}* ./data/

# write gene meta file
#zcat ${prefix}*.bed.gz | \
#cut -f 1-4 | tail -n+2 | \
#awk 'BEGIN {printf("#chr\tstart\tend\tgene\n")} {print $0}' > \
#gene.meta.txt

Rscript ${dir_GTEx}'prep.exp.R' ./data/${prefix}*.bed.gz ./data/'ex.rds'

rm ./data/${prefix}*

if ! [ -d "./script" ]; then
    mkdir "./script"
fi

if ! [ -d "./logs" ]; then
    mkdir "./logs"
fi

if ! [ -d "./plots" ]; then
    mkdir "./plots"
fi

done
