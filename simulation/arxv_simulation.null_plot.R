rm(list = ls())
parameter = commandArgs(trailingOnly = T)
file_dat_null = "./new_Sigma/simulation.null.lambda0.1.K101.rds"
file_qqplot_null = './new_Sigma/qqplot_null.png'
if_plot_qq = TRUE
if_plot_hist = FALSE

library(tidyverse)
library(gridExtra)
library("ggsci")

p_null_all = readRDS(file_dat_null)

if(if_plot_qq %in% c("True", "TRUE", "T")){
  ### qqplot with three methods on
  p_null_all = lapply(p_null_all, function(x)  -(log10(sort(x))))
  p_null_all = data.frame(p_null_all)
  
  expected <- c(1:length(p_null_all$p.null.PCO))
  lexp <- -(log10(expected / (length(expected+1))))
  
  df = cbind(p_null_all, data.frame("lexp" = lexp))
  # for small test
  #ind_sample = c(sample(1:length(p_null_all[[1]]), 1e+3), 1:100, (10^7-100):10^7)
  #df = df[sort(ind_sample), ]
  df = pivot_longer(df, !lexp, names_to = "Method", values_to = "lobs")
  
  fig <- ggplot(df, aes(x = lexp, y = lobs, color = factor(Method,
                                                           levels = c("p.null.PCO", "p.null.PC1", "p.null.minp"),
                                                           labels = c("PCO", "PC1", "minp"))
                        )) +
    geom_point(alpha = 0.7, size = 0.6) +
    geom_smooth(alpha = 0.5, size = 0.7, se = FALSE) +
    geom_abline(slope = 1, intercept = 0, color="gray") +
    #scale_color_startrek() +
    #scale_color_viridis(option = "D", discrete = TRUE) +
    scale_color_uchicago() +
    labs(x = "Expected (-logP)", y = "Observed (-logP)", color = "Method") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "right") +
    theme_classic() +
    coord_fixed(ratio = 7/8.5)
  ggsave(file_qqplot_null, fig, width = 7, height = 7)
  
}else if(if_plot_hist %in% c("True", "TRUE", "T")){
  ### QQ plot with histogram
  
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
    res[[2]] = ggplot(df, aes(yy)) + geom_histogram(breaks = seq(0, 1, 0.01)) +
      theme(text=element_text(size=6)) +
      labs(title = title, x = "Observed (P)")
    
    return(res)
  }
  
  fig1 = qqplot.hist(p_null_all$'p.null.PCO', 'p.null.PCO')
  fig2 = qqplot.hist(p_null_all$'p.null.PC1', 'p.null.PC1')
  fig3 = qqplot.hist(p_null_all$'p.null.minp', 'p.null.minp')
  fig.all = cbind(fig1, fig2, fig3)
  ggsave(file_qqplot_null,
         marrangeGrob(fig.all, ncol=length(fig.all), nrow=2, top = NULL),
         width = 4*length(fig.all), height = 4*2)
}else{cat("No plots generated.")}

