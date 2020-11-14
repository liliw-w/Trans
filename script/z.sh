plink_prefix_path=$1
expression_bed=$2
covariates_file=$3
prefix=$4
dir_script=$5
out=$6

module load R/3.6.1
#source activate /home/liliw1/.conda/envs/venv

nGene=$(($(zcat ${expression_bed} | wc -l)-1))
nBatch=$(( nGene / 100 ))
nLeft=$(( nGene % 100 ))

if (( nBatch >= 1 )); then
  if ((nLeft != 0)); then ((nBatch++)); fi
  
  for i in `seq 1 $nBatch`; do
    prefix_tmp="$prefix.$i"
    expression_bed_tmp="$expression_bed.$i.bed.gz"
    
    ~/.conda/envs/venv/bin/python3 -m tensorqtl \
    ${plink_prefix_path} ${expression_bed_tmp} ${prefix_tmp} \
    --covariates ${covariates_file} \
    --mode trans \
    --pval_threshold 1 \
    --maf_threshold 0 \
    --output_text
    
    zcat $prefix_tmp'.trans_qtl_pairs.txt.gz' | \
    tail -n+2 | \
    awk '{printf "%s\t%s\t%.6f\n", $1,$2,$4/$5}' \
    > $prefix_tmp'.trans_qtl_pairs_z.txt'
    rm $prefix_tmp'.trans_qtl_pairs.txt.gz'
    rm $prefix_tmp'.tensorQTL.trans.log'
  done
  
  cat $prefix'.'*'.trans_qtl_pairs_z.txt' | awk 'BEGIN{print "snp\tgene\tzscore"}{print $0}' > \
  $prefix'.trans_qtl_pairs_z.txt'
  rm $prefix'.'*'.trans_qtl_pairs_z.txt'
  
else
  ~/.conda/envs/venv/bin/python3 -m tensorqtl \
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
fi

# convert z.txt to z matrix
Rscript --no-save --no-restore $dir_script'make.zmat.R' $prefix".trans_qtl_pairs_z.txt" $out
rm $prefix'.trans_qtl_pairs_z.txt'

#conda deactivate

