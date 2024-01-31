# Trans-PCO Demo

## Demo code dir

See folder `/project2/xuanyao/llw/trans_pco_demo`.

## Main steps

- Below are main steps to run the pipeline. All scripts should be ready to run directly.
    - All inputs and outputs are provided for use.
    - Only main scripts are included in the demo. If you’d like scripts for any other analysis in the paper, feel free to let me know.
- I have all inputs and outputs ready. If you want to run on your own data, you can refer to the format of the input files and prepare your own data accordingly.
- It’s not necessary for you to run scripts exactly like what I did. But taking a look at the scripts, input, output files, can help understand the whole pipeline, so that you can modify the scripts for your own data/usage.
- The aim of this demo is to allow you run the code immediately and get a general sense of what inputs trans-PCO requires to run and what output of trans-PCO look like.
- For simplicity, the demo is to calculate p-values for only SNPs on chromosome 20 and gene module 150.

| Step | Goal | Script | Input (demo) | Output (demo) | Note |
| --- | --- | --- | --- | --- | --- |
| 1 | Regress out covariates | script/regress.out.covariates.R | data/ex.rds <br> data/gencode.v19.annotation.table.txt <br> data/covariates.txt | result/gene.meta.txt <br> result/ex.var.regressed.rds |  |
| 2 | Define gene modules | script/coexp.module.R | result/ex.var.regressed.rds | result/coexp.module.rds <br> result/Nmodule.txt | You may want to define your own gene modules other than using WGCNA (Require R package: WGCNA). You could look at the format of coexp.module.rds and make your own modules. |
| 3.1 | Prepare bed file of a module to run tensorqtl | script/prep.bed.R | data/ex.rds <br> result/gene.meta.txt <br> result/coexp.module.rds | result/expression.module150.bed.gz |  |
| 3.2 | Calculate z-scores using tensorqtl | script/z.sh | data/chr20_QCed <br> result/expression.module150.bed.gz <br> data/covariates.txt | z/z.module150.chr20.txt.gz | Require package: tensorqtl  <br> You don’t have run tensorqtl if you already have your own z-scores. You could check z/z.module150.chr20.txt.gz for z-score file format. |
| 4* | Calculate pvalue by PCO (trans genes) | script/p.R | result/ex.var.regressed.rds <br> result/gene.meta.txt <br> result/coexp.module.rds  <br> z/z.module150.chr20.txt.gz | p/p.module150.chr20.rds | This is the main test script using PCO.  |
| 5 | Claim significance |  |  |  | I used permutations to correct for the pvalues and define significance under FDR<0.05. But permutations can take a long time. For simplicity, you could first try using the simple way - Bonferroni correction, i.e. p*#number_of_test for each gene module. |
|  |  |  |  |  |  |


