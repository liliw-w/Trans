require(mvtnorm)
require(ggplot2)
library(ggpubr)

qqplot.hist <- function(input, title){
  observed <- sort(input)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected+1))))
  df = data.frame(x = lexp, y = lobs, yy = observed)

  res = list()
  res[[1]] = ggplot(df, aes(x=x, y=y)) + geom_point(size = 0.2) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    theme(text=element_text(size=6)) +
    labs(title = title,
         x = "Expected (-logP)", y = "Observed (-logP)") + xlim(0, 7) + ylim(0, 7)
  res[[2]] = ggplot(df, aes(yy)) + geom_histogram(breaks = seq(0, 1, 0.05)) +
    theme(text=element_text(size=6)) +
    labs(title = title, x = "Observed (P)")

  return(res)
}


params1 = '../DGN_WGCNA/script/'
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))

N = 10^6

Sigma = as.matrix(readRDS('Sigma.Muscle_Skeletal.module15.chr7.rds'))
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors


res.null = matrix(rep(NA, 8), ncol = 4,
             dimnames = list(c('alpha=0.05', 'alpha=5e-8'), c('alpha', "PC1", "PCO", "MinP")))
res.null[, 'alpha'] = c(0.05, 5e-8)
p.null.all  = list()

# z null
Z.null = rmvnorm(N, rep(0, K), Sigma)

# PC1
PC1 = Z.null %*% eigen.vec[, 1]
p.null.PC1 = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
p.null.all$'p.null.PC1' = p.null.PC1
res.null['alpha=0.05', "PC1"] = sum(p.null.PC1<0.05)/N
res.null['alpha=5e-8', "PC1"] = sum(p.null.PC1<5*10^(-8))/N
cat("PC1 done.")

# PCO
SigmaO = ModifiedSigmaOEstimate(Sigma)
p.null.PCO = ModifiedPCOMerged(Z.mat=Z.null, Sigma=Sigma, SigmaO=SigmaO)
p.null.all$'p.null.PCO' = p.null.PCO
res.null['alpha=0.05', "PCO"] = sum(p.null.PCO<0.05)/N
res.null['alpha=5e-8', "PCO"] = sum(p.null.PCO<5*10^(-8))/N
cat("PCO done.")

# MinP
#p.null.MinP = apply(Z.null, 1, function(x){
#  upper = qnorm(1 - min(pnorm(-abs(x))));
#  1 - pmvnorm(lower = -rep(upper, K), upper = rep(upper, K), mean = rep(0, K), sigma = Sigma)
#})
#p.null.all$'p.null.MinP' = p.null.MinP
#res.null['alpha=0.05', "MinP"] = sum(p.null.MinP<0.05)/N
#res.null['alpha=5e-8', "MinP"] = sum(p.null.MinP<5*10^(-8))/N
#cat("MinP done.")


# plot p.null to see if it's well calibrated, i.e. uniformly distributed
fig1 = qqplot.hist(p.null.PC1, 'p.null.PC1')
fig2 = qqplot.hist(p.null.PCO, 'p.null.PCO')
#fig3 = ggplot(as.data.frame(p.null.all), aes(x=-log10(p.null.MinP))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20))
fig.all = cbind(fig1, fig2)
ggsave('qqplot.null.png',
       marrangeGrob(fig.all, ncol=length(fig.all)/2, nrow=2, top = NULL),
       width = 4, height = 4)


print(res.null)
save(res.null, p.null.all, Z.null, file = "simulation.null.RData")


# Z alternative

require(mvtnorm)

params1 = '~/xuanyao_llw/GTEx_v8/Muscle_Skeletal/script/'
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))

N.sample = 10^4; N.sim = 10^4

Sigma = as.matrix(readRDS('Sigma.DGN.module13_chr3.100.rds'))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors

# various causal gene percentage under fixed total variance
{
  N = 500; var.b.fix = 0.2
  caus.seq = c(1, 10, 20, 30, 50, 70, 100)
  models = paste0("caus.fix=", caus.seq)

  res.alt = NULL
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    B[, models] = sapply(caus.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b.fix/x)), rep(0, K-x))) ) * sqrt(N)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% solve(Sigma) %*% B[, model])
      p.alt.Oracle = 1 - pnorm(T.oracle/as.numeric(sqrt(t(B)[model, ] %*% solve(Sigma) %*% B[, model])))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.Oracle<0.05)/N.sample))

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p.alt.PC1 = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PC1<0.05)/N.sample))

      # PCO
      method.tmp = "PCO"
      p.alt.PCO = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PCO<0.05)/N.sample))

      method.tmp = "minp"
      p.alt.minp = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp<0.05)/N.sample))

      method.tmp = "minp2"
      p.alt.minp2 = apply(Z.alt, 1, function(x) 1-(1-min(1-pchisq(x^2, 1)))^length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp2<0.05)/N.sample))

    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      colnames(res.alt) = c("model", "method", "power")
      saveRDS(res.alt, "simulation.alt.caus_fix.rds")}
  }

}


# various causal gene percentage
{
  N = 500; var.b = 0.02
  caus.seq = c(1, 10, 30, 50, 70, 100)
  models = paste0("caus=", caus.seq)

  res.alt = NULL
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    B[, models] = sapply(caus.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b)), rep(0, K-x))) ) * sqrt(N)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% solve(Sigma) %*% B[, model])
      p.alt.Oracle = 1 - pnorm(T.oracle/as.numeric(sqrt(t(B)[model, ] %*% solve(Sigma) %*% B[, model])))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.Oracle<0.05)/N.sample))

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p.alt.PC1 = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PC1<0.05)/N.sample))

      # PCO
      method.tmp = "PCO"
      p.alt.PCO = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PCO<0.05)/N.sample))

      method.tmp = "minp"
      p.alt.minp = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp<0.05)/N.sample))

      method.tmp = "minp2"
      p.alt.minp2 = apply(Z.alt, 1, function(x) 1-(1-min(1-pchisq(x^2, 1)))^length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp2<0.05)/N.sample))

    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      colnames(res.alt) = c("model", "method", "power")
      saveRDS(res.alt, "simulation.alt.caus.rds")}
  }
}


# various sample size N
{
  # Different true effects
  var.b = 0.02
  N.seq = c(200, 400, 600, 800, 1000)
  models = paste0("N=", N.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")
  #Methods = c("Oracle", "PC1", "PCO", "MinP")


  res.alt = NULL
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    beta.n = as.matrix(rnorm(K, sd = sqrt(var.b)))
    B[, models] = beta.n %*% sqrt(N.seq)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% solve(Sigma) %*% B[, model])
      p.alt.Oracle = 1 - pnorm(T.oracle/as.numeric(sqrt(t(B)[model, ] %*% solve(Sigma) %*% B[, model])))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.Oracle<0.05)/N.sample))

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p.alt.PC1 = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PC1<0.05)/N.sample))

      # PCO
      method.tmp = "PCO"
      p.alt.PCO = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PCO<0.05)/N.sample))

      method.tmp = "minp"
      p.alt.minp = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp<0.05)/N.sample))

      method.tmp = "minp2"
      p.alt.minp2 = apply(Z.alt, 1, function(x) 1-(1-min(1-pchisq(x^2, 1)))^length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp2<0.05)/N.sample))

      # MinP
      #p.alt.MinP = apply(Z.alt, 1, function(x){
      #  upper = qnorm(1 - min(pnorm(-abs(x))));
      #  1 - pmvnorm(lower = -rep(upper, K), upper = rep(upper, K), mean = rep(0, K), sigma = Sigma)
      #})
      #tmp = paste0('p.alt.MinP.', model); p.alt.all[[tmp]] = p.alt.MinP
      #res.alt[model, "MinP"] = sum(p.alt.MinP<0.05)/N.sample
      #cat("MinP done.")

    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      colnames(res.alt) = c("model", "method", "power")
      saveRDS(res.alt, "simulation.alt.N.rds")}
  }

}

# various beta variance
{
  # Different true effects
  N = 500
  var.seq = c(0.005, 0.01, 0.05, 0.1, 0.2)
  models = paste0("var=", var.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")
  #Methods = c("Oracle", "PC1", "PCO", "MinP")

  res.alt = NULL
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    B[, models] = sapply(var.seq, function(x) as.matrix(rnorm(K, sd = sqrt(x)))) * sqrt(N)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% solve(Sigma) %*% B[, model])
      p.alt.Oracle = 1 - pnorm(T.oracle/as.numeric(sqrt(t(B)[model, ] %*% solve(Sigma) %*% B[, model])))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.Oracle<0.05)/N.sample))

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p.alt.PC1 = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PC1<0.05)/N.sample))

      # PCO
      method.tmp = "PCO"
      p.alt.PCO = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.PCO<0.05)/N.sample))

      method.tmp = "minp"
      p.alt.minp = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp<0.05)/N.sample))

      method.tmp = "minp2"
      p.alt.minp2 = apply(Z.alt, 1, function(x) 1-(1-min(1-pchisq(x^2, 1)))^length(x) )
      res.alt = rbind(res.alt, c(model, method.tmp, sum(p.alt.minp2<0.05)/N.sample))


      # MinP
      #p.alt.MinP = apply(Z.alt, 1, function(x){
      #  upper = qnorm(1 - min(pnorm(-abs(x))));
      #  1 - pmvnorm(lower = -rep(upper, K), upper = rep(upper, K), mean = rep(0, K), sigma = Sigma)
      #})
      #tmp = paste0('p.alt.MinP.', model); p.alt.all[[tmp]] = p.alt.MinP
      #res.alt[model, "MinP"] = sum(p.alt.MinP<0.05)/N.sample
      #cat("MinP done.")
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      colnames(res.alt) = c("model", "method", "power")
      saveRDS(res.alt, "simulation.alt.var.rds")}
  }

}


modify <- function(x){
  y = as.data.frame(readRDS(x))
  y$power = as.numeric(levels(y$power))[y$power]
  y$model = factor(y$model, levels = unique(y$model))
  y$method = factor(y$method, levels = unique(y$method))
  return(y)
}

res.alt.N = modify('simulation.alt.N.rds')
res.alt.N2 = modify('simulation.alt.N2.rds')

res.alt.var = modify('simulation.alt.var.rds')
res.alt.var2 = modify('simulation.alt.var2.rds')

res.alt.caus = modify('simulation.alt.caus.rds')

res.alt.caus_fix = modify('simulation.alt.caus_fix.rds')

#colnames(res.alt.var) = c("model", "method", "power")

# plots

fig1 = ggplot(res.alt.var, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw()
fig2 = ggplot(res.alt.var2, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw()
fig3 = ggplot(res.alt.caus, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw()
fig4 = ggplot(res.alt.caus_fix, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw()

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

df.alt.N <- summarySE(res.alt.N, measurevar="power", groupvars=c("model", "method"))
df.alt.var <- summarySE(res.alt.var, measurevar="power", groupvars=c("model", "method"))

df.alt.caus_fix <- summarySE(res.alt.caus_fix, measurevar="power", groupvars=c("model", "method"))

fig3 = ggplot(df.alt.N, aes(x=model, y=power, group=method, color=method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=power-sd, ymax=power+sd), width=.3,
                position=position_dodge(0.05)) + theme_bw()
fig4 = ggplot(df.alt.var, aes(x=model, y=power, group=method, color=method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=power-sd, ymax=power+sd), width=.3,
                position=position_dodge(0.05)) + theme_bw()
fig5 = ggplot(df.alt.N, aes(x=model, y=power, group=method, color=method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.7,
                position=position_dodge(0.05)) + theme_bw()
fig6 = ggplot(df.alt.var, aes(x=model, y=power, group=method, color=method)) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.7,
                position=position_dodge(0.05)) + theme_bw()

ggplot(df.alt.caus_fix, aes(x=model, y=power, group=method, fill=method)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.7,
                position=position_dodge()) +
  theme_bw() +
  scale_fill_manual(values=wes_palette(name="Royal2", 5))


ggsave("power.png", ggarrange(fig1, fig2, fig3, fig4, fig5, fig6,
                              ncol = 2, nrow = 3, common.legend=T, labels = LETTERS[1:6]),
       width = 10, height = 5)


B[, "u1"] = 10*eigen.vec[, 1]
B[, "uPCO"] = 4*eigen.vec[, sum(lambdas>1)]
B[, "uK"] = 1.5*eigen.vec[, K]
