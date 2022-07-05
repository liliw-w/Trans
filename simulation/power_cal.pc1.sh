module load R/3.6.1
cd /project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/


Rscript --no-restore --no-save /home/liliw1/Trans/simulation/power_cal.R \
simulation.alt.pc1.rds \
simulation.null.lambda0.1.K101.rds \
power.pc1.rds \
0.1 \
FALSE

