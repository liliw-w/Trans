###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Run PCO on the simulated z's ###########
########### Using the nullz-based sigma ###########
###############################################################
module <- as.numeric(snakemake@params[['module']])


# I/O & paras -----
PCO.script <- '~/xuanyao_llw/DGN_no_filter_on_mappability/script/'
file_Sigma <- paste0("/project2/xuanyao/llw/eQTLGen_est_Sigma/inflation/Sigma_DGN_module", module, ".rds")
file_Z <- paste0('/project2/xuanyao/llw/eQTLGen_est_Sigma/inflation/Z_', basename(file_Sigma))

file_Sigma_nullz_seq <- list.files(path = "inflation",
                                   pattern = paste0("^Sigma_nullz\\d+\\_", basename(file_Sigma)),
                                   full.names = TRUE)

## output -----


# library PCO script -----
source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))


# read files -----
Z <- readRDS(file_Z)
rownames(Z) <- paste0('snp', 1:nrow(Z))


# run PCO -----
for(file_Sigma_nullz in file_Sigma_nullz_seq){
  Sigma_nullz <- as.matrix(readRDS(file_Sigma_nullz))
  SigmaO_nullz <- ModifiedSigmaOEstimate(Sigma_nullz)
  
  p_nullz <- as.numeric(ModifiedPCOMerged(Z.mat=Z, Sigma=Sigma_nullz, SigmaO=SigmaO_nullz))
  
  ## write out p of nullz -----
  file_p_nullz <- paste0('inflation/p_', basename(file_Sigma_nullz))
  saveRDS(p_nullz, file_p_nullz)
  
  cat(file_Sigma_nullz, "case is done. \n\n")
}


###############################################################
############### for snakemake pipeline use ###############
###############################################################
file_snm <- paste0('inflation/p_Sigma_nullz_', basename(file_Sigma))
saveRDS(NULL, file_snm)

