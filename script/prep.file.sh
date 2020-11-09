wdir='/scratch/midway2/liliw1/TCGA'
dir_data='/project2/xuanyao/llw/breastcancerTCGA/txt_file/'

cd $wdir

# regress out covariates
Rscript --no-save --no-restore ./script/regress.out.covariates.R $dir_data

# remove low quality genes
Rscript --no-save --no-restore ./script/low.quality.genes.R $dir_data

# coexpression module
Rscript --no-save --no-restore ./script/coexp.module.R


# sample name file
cat $dir_data'chr22_QCed.fam' | cut -d" " -f1 > $dir_data'sample.name.txt'

# prepare genotype file, snp name is chr:pos

