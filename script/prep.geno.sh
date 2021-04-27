#!/bin/bash 

DIR='/project2/xuanyao/data/GTEx_V8/genotype'
geno_prefix='split_by_chrom/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased'

# split WGS by chromosomes. script from Yang's lab
#from subprocess import call
#chroms = ["chr"+str(x) for x in range(1, 23)] + ["chrX", "chrY", "chrM"]
#print(chroms)
#vcf = "GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
#for c in chroms:
#    outf = "GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.%s.vcf" % c
#    call("bcftools view -r %s %s > %s" % (c, vcf, outf), shell=True) 


# filter variants; convert vcf files to bfiles; format snp's names as chr:pos
#module load plink

for chr in `seq 22`
do
plink --vcf $DIR/${geno_prefix}.chr${chr}.vcf \
--geno 0.01 \
--maf 0.05 \
--hwe 1e-6 \
--make-bed \
--silent \
--out $DIR/chr${chr}_QCed

cat $DIR/chr${chr}_QCed.bim | \
awk '{key=sprintf("%s:%s", $1, $4); printf("chr%s\t%s\t%s\t%s\t%s\t%s\t\n", $1,key,$3,$4,$5,$6)}' > $DIR/chr${chr}_QCed_dup.bim
rm $DIR/chr${chr}_QCed.bim
mv $DIR/chr${chr}_QCed_dup.bim $DIR/chr${chr}_QCed.bim

echo chr$chr
done



