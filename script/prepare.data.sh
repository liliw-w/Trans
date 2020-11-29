CHR=$1
DIR=$2
sample_DIR=$3
~/src/plink \
--vcf ${DIR}/$CHR.vcf.gz \
--double-id --biallelic-only strict --geno 0.01 --maf 0.05 --hwe 1e-6 --keep ${sample_DIR}plink_samples.txt --indiv-sort f ${sample_DIR}plink_samples.txt --make-bed --out ${DIR}chr${CHR}_QCed

# convert bfile dataset (.bed) to text dataset (.ped)
~/src/plink --bfile ${DIR}chr${CHR}_QCed --recode12 --out ${DIR}chr${CHR}_QCed

plink2 --make-bed --output-chr chrM --bfile {} --maf 0.05 --geno 0.01 --out {}.format(plink_prefix_path, plink_filtered_path)
plink --bfile {} --maf 0.05 --geno 0.01 --make-bed --output-chr chrM --keep-allele-order --out {}'.format(plink_prefix_path, plink_filtered_path)

module load plink

# split WGS by chromosomes


# filter variants and convert to plink bfiles
plink \

