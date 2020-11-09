module_id=1
wdir='/scratch/midway2/liliw1/TCGA/'
dir_data='/project2/xuanyao/llw/breastcancerTCGA/txt_file/'

data_type='obs'
p_method='TruncPCO'

cd $wdir

module load R/3.6.1

if [ "$data_type"="obs" ]
then
  type_file=''
elif [ "$data_type"="null" ]
then
  type_file='null.'
else
  echo 'Use a correct data_type.'
fi


echo 'This is z computation by tensorqtl for module'$module_id'. Use '$data_type
source activate /home/liliw1/.conda/envs/venv
for chr in `seq 22`
do
plink_prefix_path=$dir_data'chr'$chr'_QCed'
expression_bed=$dir_data'expression.'$type_file'module'$module_id'.bed.gz'
covariates_file=$dir_data'covariates.'$type_file'txt'
prefix='./z/module'$module_id'.chr'$chr

python3 -m tensorqtl \
${plink_prefix_path} ${expression_bed} ${prefix} \
--covariates ${covariates_file} \
--mode trans \
--pval_threshold 1 \
--maf_threshold 0 \
--output_text

zcat $prefix'.trans_qtl_pairs.txt.gz' | \
tail -n+2 | \
awk 'BEGIN{print "snp\tgene\tzscore"}{printf "%s\t%s\t%.6f\n", $1,$2,$4/$5}' \
> $prefix'.trans_qtl_pairs_z.txt'
rm $prefix'.trans_qtl_pairs.txt.gz'
rm $prefix'.tensorQTL.trans.log'

# convert z.txt to z matrix
Rscript --no-save --no-restore ./script/make.zmat.R $data_type $module_id $chr
rm $prefix'.trans_qtl_pairs_z.txt'

done
conda deactivate
echo 'z for module'$module_id' is done. Use '$data_type

