require(mvtnorm)
require(ggplot2)
require(gridExtra)

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


params1 = './script/'
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))

Sigma = as.matrix(readRDS('Sigma.DGN.module13_chr3.100.rds'))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors

N.sample = 10^7
p.null.all  = list()
res.null = NULL

# z null
Z.null = rmvnorm(N.sample, rep(0, K), Sigma)


# PCO
method.tmp = "PCO"
p.null = as.numeric(ModifiedPCOMerged(Z.mat=Z.null, Sigma=Sigma, SigmaO=SigmaO))

p.null.all$'p.null.PCO' = p.null
res.null = rbind(res.null, c('alpha=0.05', method.tmp, sum(p.null<0.05)/N.sample))
res.null = rbind(res.null, c('alpha=1e-4', method.tmp, sum(p.null<10^(-4))/N.sample))
res.null = rbind(res.null, c('alpha=1e-5', method.tmp, sum(p.null<10^(-5))/N.sample))

cat("PCO done.")

print(res.null)

save(res.null, p.null.all, file = "simulation.null.lambda1.RData")


# PC1
method.tmp = "PC1"
PC1 = Z.null %*% eigen.vec[, 1]
p.null = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

p.null.all$'p.null.PC1' = p.null
res.null = rbind(res.null, c('alpha=0.05', method.tmp, sum(p.null<0.05)/N.sample))
res.null = rbind(res.null, c('alpha=1e-4', method.tmp, sum(p.null<10^(-4))/N.sample))
res.null = rbind(res.null, c('alpha=1e-5', method.tmp, sum(p.null<10^(-5))/N.sample))

cat("PC1 done.")


# univariate minp
method.tmp = "minp"
p.null = apply(Z.null, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

p.null.all$'p.null.minp' = p.null
res.null = rbind(res.null, c('alpha=0.05', method.tmp, sum(p.null<0.05)/N.sample))
res.null = rbind(res.null, c('alpha=1e-4', method.tmp, sum(p.null<10^(-4))/N.sample))
res.null = rbind(res.null, c('alpha=1e-5', method.tmp, sum(p.null<10^(-5))/N.sample))

cat("minp done.")

print(res.null)

save(res.null, p.null.all, file = "simulation.null.lambda1.RData")

# plot p.null to see if it's well calibrated, i.e. uniformly distributed
fig1 = qqplot.hist(p.null.all$'p.null.PC1', 'p.null.PC1')
fig2 = qqplot.hist(p.null.all$'p.null.minp', 'p.null.minp')
fig3 = qqplot.hist(p.null.all$'p.null.PCO', 'p.null.PCO')
fig.all = cbind(fig1, fig2, fig3)
ggsave('qqplot.null.png',
       marrangeGrob(fig.all, ncol=length(fig.all)/2, nrow=2, top = NULL),
       width = length(fig.all), height = 4)



# plot
{
  require(ggplot2)
  library(ggpubr)

  modify <- function(x){
    y = as.data.frame(readRDS(x))
    y$power = as.numeric(levels(y$power))[y$power]
    y$model = factor(y$model, levels = unique(y$model))
    y$method = factor(y$method, levels = unique(y$method))
    return(y)
  }

  res.alt.N = modify('simulation.alt.N.lambda0.1.K105.rds')
  #res.alt.N2 = modify('simulation.alt.N2.rds')

  #res.alt.var = modify('simulation.alt.N.05.rds')
  #res.alt.var2 = modify('simulation.alt.var2.rds')

  res.alt.caus = modify('simulation.alt.caus.lambda0.1.K105.rds')
  res.alt.caus_fix = modify('simulation.alt.caus_fix.lambda0.1.K105.rds')

  #colnames(res.alt.var) = c("model", "method", "power")


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
  #df.alt.N2 <- summarySE(res.alt.N2, measurevar="power", groupvars=c("model", "method"))
  #df.alt.var <- summarySE(res.alt.var, measurevar="power", groupvars=c("model", "method"))
  #df.alt.var2 <- summarySE(res.alt.var2, measurevar="power", groupvars=c("model", "method"))
  df.alt.caus <- summarySE(res.alt.caus, measurevar="power", groupvars=c("model", "method"))
  df.alt.caus_fix <- summarySE(res.alt.caus_fix, measurevar="power", groupvars=c("model", "method"))

  # plots

  fig1 = list(ggplot(res.alt.N, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw(),
              ggplot(data = df.alt.N, aes(x=model, y=power, group=method, color=method)) +
                geom_line() +
                geom_point() +
                geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                              position=position_dodge(0)) +
                coord_cartesian(ylim=c(0, 1)) + theme_bw()
  )

  fig2 = list(ggplot(res.alt.var, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw(),
              ggplot(df.alt.var, aes(x=model, y=power, group=method, color=method)) +
                geom_line() +
                geom_point() +
                geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                              position=position_dodge(0)) +
                coord_cartesian(ylim=c(0, 1)) + theme_bw()
  )

  fig3 = list(ggplot(res.alt.caus, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw(),
              ggplot(df.alt.caus, aes(x=model, y=power, group=method, color=method)) +
                geom_line() +
                geom_point() +
                geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                              position=position_dodge(0)) +
                coord_cartesian(ylim=c(0, 1)) + theme_bw()
  )

  fig4 = list(ggplot(res.alt.caus_fix, aes(x=model, y=power, col = method)) + geom_boxplot() + coord_cartesian(ylim=c(0, 1)) + theme_bw(),
              ggplot(df.alt.caus_fix, aes(x=model, y=power, group=method, color=method)) +
                geom_line() +
                geom_point() +
                geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                              position=position_dodge(0)) +
                coord_cartesian(ylim=c(0, 1)) + theme_bw()
  )

  ggsave("power.lambda0.1.K105.png", ggarrange(plotlist = c(fig1, fig3, fig4),
                                ncol = 2, nrow = 3, common.legend=T,
                                labels = rep(c("N", "caus", "caus_fix"), each=2)),
         width = 10, height = 15)

}


B[, "u1"] = 10*eigen.vec[, 1]
B[, "uPCO"] = 4*eigen.vec[, sum(lambdas>1)]
B[, "uK"] = 1.5*eigen.vec[, K]
