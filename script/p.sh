module_id=1
wdir='/scratch/midway2/liliw1/TCGA/'
dir_data='/project2/xuanyao/llw/breastcancerTCGA/txt_file/'

data_type='obs'
p_method='TruncPCO'

cd $wdir

echo 'This is p computation by '$p_method' for module'$module_id'. Use '$data_type
for chr in `seq 22`
do
Rscript --no-save --no-restore ./script/p.R $data_type $module_id $chr $p_method
done
echo 'p for module'$module_id' is done. Use '$data_type
