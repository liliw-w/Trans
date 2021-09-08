rm(list = ls())
parameter = commandArgs(trailingOnly = T)
PCO.script = parameter[1] # './script_lambda0.1/'
file.Sigma = parameter[2] # 'Sigma.DGN.module13_chr3.100.rds'
file.res = parameter[3] # "simulation.null.lambda0.1.K100.rds"
is_plot = parameter[4]

require(mvtnorm)
source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))

Sigma = as.matrix(readRDS(file.Sigma))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors

N.sample = 10^7
p.null.all  = list()

# z null
Z.null = rmvnorm(N.sample, rep(0, K), Sigma)

# PCO
p.null.all$'p.null.PCO' = as.numeric(ModifiedPCOMerged(Z.mat=Z.null, Sigma=Sigma, SigmaO=SigmaO))

cat("PCO done.")
saveRDS(p.null.all, file = file.res)

# PC1
PC1 = Z.null %*% eigen.vec[, 1]
p.null.all$'p.null.PC1' = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

cat("PC1 done.")

# univariate minp
p.null.all$'p.null.minp' = apply(Z.null, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

cat("minp done.")
saveRDS(p.null.all, file = file.res)

if(is_plot %in% c("True", "TRUE", "T")){
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
           x = "Expected (-logP)", y = "Observed (-logP)")
    #res[[2]] = ggplot(df, aes(yy)) + geom_histogram(breaks = seq(0, 1, 0.05)) +
    #  theme(text=element_text(size=6)) +
    #  labs(title = title, x = "Observed (P)")

    return(res)
  }

  fig1 = qqplot.hist(p.null.all$'p.null.PCO', 'p.null.PCO')
  fig2 = qqplot.hist(p.null.all$'p.null.PC1', 'p.null.PC1')
  fig3 = qqplot.hist(p.null.all$'p.null.minp', 'p.null.minp')
  fig.all = cbind(fig1, fig2, fig3)
  ggsave('qqplot.null.png',
         marrangeGrob(fig.all, ncol=length(fig.all), nrow=1, top = NULL),
         width = length(fig.all), height = 4)
}
