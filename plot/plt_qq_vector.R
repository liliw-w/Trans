qqplot <- function(input,
                   ci_level = 0.95,
                   my_theme = "/project2/xuanyao/llw/Trans/plot/theme_my_pub.R"){
  require(tidyverse)
  source(my_theme)
  
  
  #reset obs p-values that are 0 to a fixed value, here I use the min non-zero p/10
  input[input == 0] = min(input[input != 0])/10

  # number of samples
  n = length(input)
  
  # expected
  expected = seq(1, n) / (n+1)
  lexp = -log10(expected)
  
  # order statistic of null p
  ci_l = -log10( qbeta(p = (1 - ci_level) / 2, shape1 = 1:n, shape2 = n:1) )
  ci_r = -log10( qbeta(p = (1 + ci_level) / 2, shape1 = 1:n, shape2 = n:1) )
  
  
  # obs
  observed = sort(input)
  lobs = -log10(observed)
  
  
  
  # data for plt
  df_plt = data.frame(x = lexp, y = lobs, yy = observed, ci_l = ci_l, ci_r = ci_r)
  
  return(
    ggplot(df_plt, aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = ci_l, ymax = ci_r), fill = "grey80", color="grey80") +
      geom_abline(slope = 1, intercept = 0, color = "black") +
      geom_point(aes(color = y), size = 0.3) +
      labs(x = bquote(Expected -log[10]~italic((P))), y = bquote(Observed -log[10]~italic((P)))) +
      scale_color_gradientn(colors = RColorBrewer::brewer.pal(8, "Reds")[3:8]) +
      theme_my_pub(legend.position = "none")
  )
}

