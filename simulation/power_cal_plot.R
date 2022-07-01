rm(list = ls())

### Plot
library(data.table)
library(tidyverse)
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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
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
    geom_pointrange(aes(ymin = power-ci, ymax = power+ci),
                    fatten = 1, size = 0.5,
                    position = position_dodge(width = 0.3)) +
    labs(x = xlab_name, y = "Power", color = "Method") +
    scale_colour_manual(#values = cbp2,
      breaks = c("trans-PCO", "PC1", "minp"),
      labels = c("trans-PCO" = "Trans-PCO", "minp" = "MinP", "PC1" = "PC1"),
      values = c("trans-PCO" = "#85192d", "minp" = "#e89c31", "PC1" = "#1d349a"),
      guide = guide_legend(override.aes = list(size = 0.3))
    ) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
      
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      #legend.background = element_rect(color = "black", linetype = "dashed"),
      legend.key.size= unit(0.5, "cm"),
      
      axis.line = element_line(colour="black"),
      #axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      axis.text=element_text(colour = "black", size=12),
      axis.title.y = element_text(angle=90,vjust =2, size=14),
      axis.title.x = element_text(vjust = -0.2, size=14),
      
      plot.margin=unit(c(10,5,5,5),"mm")
    ) +
    coord_cartesian(xlim = c(1.2, n_distinct(df.alt$model)-0.2))
  
  ggsave(paste0(file_out, "3.pdf"), fig, height = 3, width = 4.5)
  
  
  
  #########################################################
  #########################################################
  #########################################################
  ### 3. point range plot
  dat2 <- res.alt %>% group_by(method, model, change) %>% summarise(m = mean(power), sd = sd(power) ) %>% ungroup()
  dat2$method <- factor(dat2$method,
                          levels = c("PCO", "minp", "PC1"),
                          labels = c("trans-PCO", "minp", "PC1"))
  dat2$model = factor(dat2$model,
                        levels = levels(dat2$model),
                        labels = sapply(levels(dat2$model), function(x) strsplit(x, "=")[[1]][2] )
  )
  #dat2$change <- "Sample Size"
  dat2$change <- factor(dat2$change,
                        levels = c("Sample Size", "Causality", "Genetic Variance"),
                        labels = c("Sample Size", "Causality Proportion", "Genetic Variance"))
  
  ggplot(dat2, aes(x = model, y = m, color = method)) +
    #geom_point(data = res.alt,
    #           aes(x = factor(model), y = power, color = method),
    #           position = position_jitterdodge(), alpha = 0.1) +
    geom_pointrange(aes(ymin = m-sd, ymax = m+sd),
                    fatten = 2.5, size = 0.7,
                    position = position_dodge(width = 0.3)) +
    labs(x = "Model", y = "Power", color = "Method") + #x = xlab_name, 
    facet_wrap(vars(change), scales = "free_x") +
    scale_colour_manual(#values = cbp2,
      breaks = c("trans-PCO", "minp", "PC1"),
      labels = c("Trans-PCO", "Minp", "PC1-based"),
      values = #c("trans-PCO" = "#074e67", "PC1" = "#67074e", "minp" = "#dd9933")
        c("trans-PCO" = "#85192d", "minp" = "#e89c31", "PC1" = "#1d349a"),
      guide = guide_legend(override.aes = list(size = 0.3))
    ) +
    theme_classic(base_size = 16) +
    theme(
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
          
          legend.position = "bottom", # c(0.25, 0.8),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          legend.background = element_rect(color = "black", linetype = "dashed"),
          legend.key.size= unit(0.5, "cm"),
          
          axis.line = element_line(colour="black"),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_blank(),
          
          
          plot.margin=unit(c(10,5,5,5),"mm"),
          axis.text=element_text(colour = "black", size=14),
          axis.title.y = element_text(angle=90,vjust =2, size=16),
          axis.title.x = element_text(vjust = -0.2, size=16) )
  
  
  ggsave("updated_power_plot.pdf", height = 4, width = 9)
  
  dat2[dat2$change=="N", "change"] <- "Sample Size"
  
  
  
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
