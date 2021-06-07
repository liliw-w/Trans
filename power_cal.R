rm(list = ls())
parameter = commandArgs(trailingOnly = T)
file_p_alt = parameter[1] # "simulation.alt.N.lambda0.1.K100.rds"
file_p_null = parameter[2] # "simulation.null.lambda0.1.K100.rds"
file_out = parameter[3] # "power.N.lambda0.1.K100.rds"
is_plot = parameter[4]

require(data.table)
p_alt_all = readRDS(file_p_alt)
p_null_all = readRDS(file_p_null)

N.sample = dim(p_alt_all)[1]
C = length(p_null_all[[1]])/N.sample
res = NULL
for (method.tmp in c("PCO", "PC1", "minp")) {
  p.null = data.table(paste0("null", 1:length(p_null_all[[paste0("p.null.", method.tmp)]])), p_null_all[[paste0("p.null.", method.tmp)]])
  res[[method.tmp]] = apply(p_alt_all[ , method.tmp, , ], c(2, 3), function(x) {
    p.obs = data.table(paste0("obs", 1:N.sample), as.numeric(x) )
    p.obs.rank = frank(p.obs, V2)
    names(p.obs.rank) = p.obs$V1

    all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
    names(all.rank) = c(p.obs$V1, p.null$V1)
    q = pmin((all.rank[names(p.obs.rank)]/p.obs.rank-1)/C, 1)

    sum(q < 0.05)/N.sample
  } )
  cat("Method ", method.tmp, "is done.", '\n')
  saveRDS(res, file_out)
}

### Plot
if(is_plot %in% c("True", "TRUE", "T")){
  require(ggplot2)
  require(ggpubr)
  modify <- function(x){
    power_emp = readRDS(x)
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
               xlab("Sample size") +
               ylab("Power") +
               labs(color = "Method") +
               geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                             position=position_dodge(0)) +
               coord_cartesian(ylim=c(0, 1)) + theme_bw() + scale_color_manual(values=c("red3", "blue2", "springgreen4"))
  )
  ggsave(paste0(file_out, ".png"), ggarrange(plotlist = c(fig),
                                             ncol = 2, nrow = 1, common.legend=T,
                                             labels = "A"),
         height = 4.5, width = 7)
}
