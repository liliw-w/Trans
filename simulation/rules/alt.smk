rule sim_alt_low_caus_high_b:
    input:
      file_Sigma='/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/Sigma-new_DGN_module29_K101.rds'
    output:
      file_p_alt='/scratch/midway3/liliw1/paper1_sim/simulation_alt_lowCaus_highb_change'+change+'_varb{varb}_N'+N+'_K101.rds'
    params:
      change=change,
      oracle_thre=oracle_thre,
      dir_pco=dir_pco,
      caus=caus,
      N=N,
      N_sample=N_sample,
      N_sim=N_sim
    shell:
      """
      {dir_env_r}/bin/Rscript --vanilla ~/Trans/simulation/2_5_sim_alt_low_caus_high_b.R {params.change} {params.oracle_thre} {params.dir_pco} {input.file_Sigma} {output.file_p_alt} {wildcards.varb} {params.caus} {params.N} {params.N_sample} {params.N_sim}
      """


rule cal_power:
    input:
      file_p_alt='/scratch/midway3/liliw1/paper1_sim/simulation_alt_lowCaus_highb_change'+change+'_varb{varb}_N'+N+'_K101.rds',
      file_dat_null="/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/simulation.null.lambda0.1.K101.rds"
    output:
      file_power='/scratch/midway3/liliw1/paper1_sim/power_lowCaus_highb_change'+change+'_varb{varb}_N'+N+'_K101.rds'
    params:
      fdr_level=fdr_level
    shell:
      """
      {dir_env_r}/bin/Rscript --vanilla ~/Trans/simulation/3_cal_power.R {input.file_p_alt} {input.file_dat_null} {params.fdr_level} {output.file_power}
      """

rule plt_power:
    input:
      file_power='/scratch/midway3/liliw1/paper1_sim/power_lowCaus_highb_change'+change+'_varb{varb}_N'+N+'_K101.rds'
    output:
      file_fig_power='/scratch/midway3/liliw1/paper1_sim/plt_power_lowCaus_highb_change'+change+'_varb{varb}_N'+N+'_K101.pdf'
    shell:
      """
      {dir_env_r}/bin/Rscript --vanilla ~/Trans/simulation/4_4_plt_power_low_caus_high_b.R {input.file_power} {output.file_fig_power}
      """

