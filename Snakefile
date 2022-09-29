configfile: "config.yaml"

path='result/Nmodule.txt'
with open(path) as f:
    lines = [x.rstrip() for x in f]
Nmodule=int(lines[0])

MODULE=list(range(1, Nmodule+1))
CHRS=list(range(1, config['Nchr']+1))
PERM=list(range(1, config['Nperm']+1))
fdr_thre_chr_module=config['fdr_level']

rule all:
  input:
    #expand('p/p.module{module}.chr{chr}.rds', module=MODULE, chr=CHRS)
    'postanalysis/indep.signals.chr.module.perm'+str(config['Nperm'])+'.txt'

rule prep_bed:
  input:
    file_covariates=config['dir_expression']+config['file_covariates'],
    file_gene_meta='result/'+config['file_gene_meta'],
    file_coexp_module='result/'+config['file_coexp_module'],
    file_ex=config['dir_expression']+config['file_ex']
  output:
    file_expression=expand('result/expression.module{module}.bed.gz', module=MODULE)
  params:
    data_type='obs'
  script:
    'script/'+config['script_prep_bed']

rule z:
  input:
    expression_bed='result/expression.module{module}.bed.gz',
    file_covariates=config['dir_expression']+config['file_covariates']
  output:
    file_z='z/z.module{module}.chr{chr}.txt.gz'
  params:
    plink_prefix_path=config['dir_geno']+config['geno_prefix']+'{chr}'+config['geno_suffix'],
    prefix='module{module}.chr{chr}',
    dir_script='script/'
  shell: 'bash '+'script/'+config['script_z']+' {params.plink_prefix_path} {input.expression_bed} {input.file_covariates} {params.prefix} {params.dir_script} {output.file_z}'

rule p:
  input:
    file_ex_var_regressed='result/'+config['file_ex_var_regressed'],
    file_gene_meta='result/'+config['file_gene_meta'],
    file_coexp_module='result/'+config['file_coexp_module'],
    file_z='z/z.module{module}.chr{chr}.txt.gz'
  output:
    file_p='p/p.module{module}.chr{chr}.rds'
  params:
    dir_script='script/', chr='{chr}', module='{module}'
  script:
    'script/'+config['script_p']

rule prep_bed_null:
  input:
    file_covariates=config['dir_expression']+config['file_covariates'],
    file_gene_meta='result/'+config['file_gene_meta'],
    file_coexp_module='result/'+config['file_coexp_module'],
    file_ex=config['dir_expression']+config['file_ex']
  output:
    file_expression=temp(expand('result/expression.null.module{module}.perm{perm}.bed.gz', module=MODULE, allow_missing=True)),
    file_covariates_null=temp('result/covariates.null.perm{perm}.txt')
  params:
    data_type='null'
  script:
    'script/'+config['script_prep_bed']


rule z_null:
  input:
    expression_bed='result/expression.null.module{module}.perm{perm}.bed.gz',
    file_covariates='result/covariates.null.perm{perm}.txt'
  output:
    file_z=temp('z/z.null.module{module}.chr{chr}.perm{perm}.txt.gz')
  params:
    plink_prefix_path=config['dir_geno']+config['geno_prefix']+'{chr}'+config['geno_suffix'],
    prefix='module{module}.chr{chr}.perm{perm}.null',
    dir_script='script/'
  shell: 'bash script/'+config['script_z']+' {params.plink_prefix_path} {input.expression_bed} {input.file_covariates} {params.prefix} {params.dir_script} {output.file_z}'


rule p_null:
  input:
    file_ex_var_regressed='result/'+config['file_ex_var_regressed'],
    file_gene_meta='result/'+config['file_gene_meta'],
    file_coexp_module='result/'+config['file_coexp_module'],
    file_z='z/z.null.module{module}.chr{chr}.perm{perm}.txt.gz'
  output:
    file_p='p/p.null.module{module}.chr{chr}.perm{perm}.rds'
  params:
    dir_script='script/', chr='{chr}', module='{module}'
  script:
    'script/'+config['script_p']


rule FDR_chr_module:
  input:
    file_p=expand('p/p.module{module}.chr{chr}.rds',chr=CHRS,module=MODULE),
    file_p_null=expand('p/p.null.module{module}.chr{chr}.perm{perm}.rds',chr=CHRS,module=MODULE, allow_missing=True)
  output:
    file_q='FDR/q.chr.module.perm{perm}.rds'
  script:
    'script/'+config['script_q']

rule average_perm_chr_module:
  input:
    file_q=expand('FDR/q.chr.module.perm{perm}.rds', perm=PERM)
  output:
    file_signals='FDR/signals.chr.module.perm'+str(config['Nperm'])+'.txt'
  params:
    fdr_thre=fdr_thre_chr_module
  script:
    'script/'+config['script_average_perm']


rule post_chr_module:
  input: sig='FDR/signals.chr.module.perm'+str(config['Nperm'])+'.txt'
  output: sig_uniq='postanalysis/LD.prun.in.chr.module.perm'+str(config['Nperm'])+'.txt', sig_indp='postanalysis/indep.signals.chr.module.perm'+str(config['Nperm'])+'.txt'
  params: dir_geno=config['dir_geno'], geno_prefix=config['geno_prefix'], geno_suffix=config['geno_suffix']
  shell: 'bash script/'+config['script_post']+' {input.sig} {output.sig_uniq} {output.sig_indp} {params.dir_geno} {params.geno_prefix} {params.geno_suffix}'
