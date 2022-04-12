rm(list = ls())

### Plot
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsci)

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

cbp2 <- c("#E69F00", "#009E73", "#0072B2", "#000000", "#56B4E9",
          "#F0E442", "#D55E00", "#CC79A7")

fig_all = list()
for (change in c('N', 'caus')) {
  if(change == "N"){
    xlab_name = expression(paste("Sample Size"))
    file_out = 'new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds'
  }else if(change == "caus"){
    xlab_name = expression(paste("Causal Proportion"))
    file_out = 'new_Sigma/power.caus.lambda0.1.varb1e-3.K101.rds'
  }else if(change == "var"){
    xlab_name = expression(paste("Genetic Variance")) 
    file_out = 'new_Sigma/power.var.lambda0.1.varb1e-3.K101.rds'
  }else{cat("Wrong parameter.")}
  
  
  res.alt = modify(file_out)
  df.alt <- summarySE(res.alt, measurevar="power", groupvars=c("model", "method"), na.rm=TRUE)
  df.alt$method <- factor(df.alt$method,
                          levels = c("PCO", "PC1", "minp"),
                          labels = c("trans-PCO", "PC1", "minp"))
  df.alt$model = factor(df.alt$model,
                        levels = levels(df.alt$model),
                        labels = sapply(levels(df.alt$model), function(x) strsplit(x, "=")[[1]][2] )
  )
  
  fig <- ggplot(data = df.alt, aes(x=model, y=power, group=method, color=method)) +
    geom_line(linetype = "solid", size = 1) +
    geom_point(shape = 20, size = 2) +
    geom_errorbar(aes(ymin=power-se, ymax=power+se), width=.3,
                  position=position_dodge(0), size = 1) +
    labs(x = xlab_name, y = "Power", color = "Method") +
    scale_colour_manual(values = cbp2) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.2, 0.8),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          legend.background = element_rect(color = "black", linetype = "dashed"),
          legend.key.size= unit(0.5, "cm"),
          axis.line = element_line(colour="black"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          axis.text=element_text(colour = "black", size=16),
          axis.title.y = element_text(angle=90,vjust =2, size=16),
          axis.title.x = element_text(vjust = -0.2, size=16) ) +
    coord_fixed(ratio = nlevels(df.alt$model), ylim=c(0, 0.9))
  
  fig
  
  ggsave(paste0(file_out, ".png"), fig, height = 5, width = 5)
  
  cat("Plot stored in", paste0(file_out, ".png"), '\n')
  fig_all[[change]] = fig
  
  if(change == 'caus'){
    ### draw zoomed part
    fig_zoomed <- fig +
      coord_fixed(ratio = nlevels(df.alt$model)/0.03, xlim = c(1.5, 2.5)) +
      xlim( levels(df.alt$model)[1:3] ) +
      ylim(c(0, 0.009))
    fig_zoomed
    ggsave(paste0(file_out, "_zoomed.png"), fig_zoomed, height = 4, width = 5)
    
  }
}


### combine two plot into figure 1
ggsave('new_Sigma/Figure1_combined.png',
       ggarrange(plotlist = fig_all,
                 ncol = length(fig_all), nrow = 1,
                 labels = LETTERS[1:length(fig_all)],
                 label.y = 0.9),
       height = 5, width = 5*length(fig_all))
