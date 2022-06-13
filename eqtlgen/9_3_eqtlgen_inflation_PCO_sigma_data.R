###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Run PCO on the simulated z's ###########
########### Using the data-based sigma ###########
###############################################################

### I/O
module <- as.numeric(snakemake@params[['module']])

PCO.script <- '~/xuanyao_llw/DGN_no_filter_on_mappability/script/'
#file_Sigma <- list.files(path = "inflation",
#                         pattern = paste0("^Sigma_DGN_module", module, "\\_K\\d*\\.rds$"),
#                         full.names = TRUE)
file_Sigma <- list.files(path = "inflation",
                         pattern = paste0("^Sigma_DGN_module", module, "\\.rds$"),
                         full.names = TRUE)
file_Z <- paste0('inflation/Z_', basename(file_Sigma))
file_p <- paste0('inflation/p_', basename(file_Sigma))


### library PCO script
source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))


### read data
Sigma <- as.matrix(readRDS(file_Sigma))
Z <- readRDS(file_Z)
rownames(Z) <- paste0('snp', 1:nrow(Z))


### run PCO
SigmaO <- ModifiedSigmaOEstimate(Sigma)
p <- as.numeric(ModifiedPCOMerged(Z.mat=Z, Sigma=Sigma, SigmaO=SigmaO))


saveRDS(p, file = file_p)

