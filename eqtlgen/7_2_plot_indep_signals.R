###############################################################
########### plot independent significant SNPs ###########
########### Do we really need independent SNPs? ###########
########### as each eQTLGen SNPs here are selected disease associated SNPs ###########
###############################################################


### refer to the similar figure of DGN signals ###


rm(list = ls())
library(data.table)
library(tidyverse)
source("~/Trans/followup/theme_my_pub.R")


file_signal <- 'postanalysis/signal.txt'


### read data
signal <- fread(file_signal)


