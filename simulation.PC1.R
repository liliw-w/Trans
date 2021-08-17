# Z alternative
# FDR based on empirical pnull
rm(list = ls())
parameter = commandArgs(trailingOnly = T)

change = 'N'
oracle.thre = 0.1
PCO.script = './script_lambda0.1/'
file.Sigma = 'Sigma.DGN.module13_chr3.100.rds'
file.res = 'simulation.alt.N.lambda0.1.varb1e-3.K100.PC1.rds'
var.b = 0.001
caus = 0.3
N = 500

require(mvtnorm)

source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))


N.sample = 10^3; N.sim = 1
#thre = 0.05/(50*1000000)

Sigma = as.matrix(readRDS(file.Sigma))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors
Sigma.trunc.inv = eigen.vec[, which(lambdas>oracle.thre)] %*% diag(1/lambdas[lambdas>oracle.thre]) %*% t(eigen.vec[, which(lambdas>oracle.thre)])


  cat(paste0("Various ", change), '\n')

  # Different true effects
  N.seq = c(200, 400, 600, 800)
  models = paste0("N=", N.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")


  res.alt = NULL
  p_alt_all = array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    #beta.n = as.matrix(rnorm(K, sd = sqrt(var.b)))
    #beta.n = rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
    beta.n = as.matrix(eigen.vec[, 1] * sqrt(var.b * K))
    B[, models] = beta.n %*% sqrt(N.seq)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] = 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

      # univariate minp
      method.tmp = "minp"
      p_alt_all[, method.tmp, model, i] = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

      # PCO
      method.tmp = "PCO"
      p_alt_all[, method.tmp, model, i] = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    cat("Simulation: ", i, "\n")
    #if(i %% 20 == 0){
    #  saveRDS(p_alt_all, file.res)}
  #}
    saveRDS(p_alt_all, file.res)
  }


###### Calculate power ##########
file_p_alt = 'simulation.alt.N.lambda0.1.varb1e-3.K100.PC1.rds'
file_p_null = "simulation.null.lambda0.1.K100.rds"
fdr_level = 0.05
file_out = paste0("power.N.lambda0.1.varb1e-3.K100.PC1.fdr", fdr_level, ".rds")
is_plot = 'True'
xlab_text = "Sample Size"

require(data.table)
p_alt_all = readRDS(file_p_alt)
p_null_all = readRDS(file_p_null)

N.sample = dim(p_alt_all)[1]
C = length(p_null_all[[1]])/N.sample
res = NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  p.null = data.table(paste0("null", 1:length(p_null_all[[paste0("p.null.", method.tmp)]])), p_null_all[[paste0("p.null.", method.tmp)]])
  res[[method.tmp]] = apply(p_alt_all[ , method.tmp, , ], c(2), function(x) {
    p.obs = data.table(paste0("obs", 1:N.sample), as.numeric(x) )
    p.obs.rank = frank(p.obs, V2)
    names(p.obs.rank) = p.obs$V1

    all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
    names(all.rank) = c(p.obs$V1, p.null$V1)
    q = pmin((all.rank[names(p.obs.rank)]/p.obs.rank-1)/C, 1)

    sum(q < fdr_level)/N.sample
  } )
  cat("Method ", method.tmp, "is done.", '\n')
  print(res)
}
print(res)
saveRDS(res, file_out)


### Plot ###
if(is_plot %in% c("True", "TRUE", "T")){
  require(ggplot2)
  require(ggpubr)
  modify <- function(x){
    power_emp = readRDS(x)
    power_emp = lapply(power_emp, function(x) as.matrix(x))
    DT = data.table("method" = rep(names(power_emp), each=nrow(power_emp[[1]])),
                    "model" = rownames(do.call(rbind, power_emp)),
                    do.call(rbind, power_emp))
    y = melt(DT, id.vars = c("method", "model"), value.name = "power")[, -"variable"]

    y$power = as.numeric(y$power)
    y$model = factor(y$model, levels = unique(y$model))
    y$method = factor(y$method, levels = unique(y$method))
    return(y)
  }
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }

  res.alt = modify(file_out)
  df.alt <- summarySE(res.alt, measurevar="power", groupvars=c("model", "method"), na.rm=TRUE)
  fig = list(ggplot(data = df.alt, aes(x=model, y=power, group=method, color=method)) +
               geom_line(linetype = "solid", size = 0.15) +
               geom_point(shape = 20) +
               xlab(xlab_text) +
               ylab("Power") +
               scale_y_continuous(labels = scales::percent) +
               coord_cartesian(ylim=c(0, 1)) +
               labs(color = "Method") +
               geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                             position=position_dodge(0)) +
               theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               scale_color_manual(values=c("red3", "blue2", "springgreen4"))
  )
  ggsave(paste0(file_out, ".png"), ggarrange(plotlist = c(fig),
                                             ncol = 1, nrow = 1, labels = "C"),
         height = 5, width = 6)
}

