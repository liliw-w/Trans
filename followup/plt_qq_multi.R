qqplot <- function(input,
                   ci_level = 0.95,
                   group_title = NULL,
                   group_order = "Type",
                   is_group_numerical_order = TRUE,
                   my_theme = "/home/liliw1/Trans/followup/theme_my_pub.R"){
  require(tidyverse)
  source(my_theme)
  
  #reset obs p-values that are 0 to a fixed value, here I use the min non-zero p/10
  input[input == 0] = min(input[input != 0])/10
  
  # number of samples
  n = nrow(input)
  
  # expected
  expected = seq(1, n) / (n+1)
  lexp = -log10(expected)
  
  # order statistic of null p
  ci_l = -log10( qbeta(p = (1 - ci_level) / 2, shape1 = 1:n, shape2 = n:1) )
  ci_r = -log10( qbeta(p = (1 + ci_level) / 2, shape1 = 1:n, shape2 = n:1) )
  
  
  # obs
  observed = apply(input, 2, sort) %>% as.data.frame()
  lobs = -log10(observed)
  
  
  # data for plt
  df_plt = cbind(data.frame(x = lexp, ci_l = ci_l, ci_r = ci_r), lobs) %>%
    pivot_longer(-c(x, ci_l, ci_r), names_to = "Type", values_to = "y")
  # set group order in plt
  group_order = if(is.null(group_order)) unique(df_plt$Type) else group_order
  df_plt$Type = factor(df_plt$Type, levels = group_order)
  
  if(is_group_numerical_order){
    ggplot(df_plt, aes(x = x, y = y, group = Type)) +
      geom_ribbon(aes(ymin = ci_l, ymax = ci_r), fill = "grey80", color="grey80") +
      geom_abline(slope = 1, intercept = 0, color = "black", size = 1) +
      geom_point(aes(color = Type, shape = Type), size = 0.5) +
      labs(x = bquote(Expected -log[10]~italic((P))),
           y = bquote(Observed -log[10]~italic((P))),
           color = group_title,
           shape = group_title) +
      scale_shape_manual(values = 1:nlevels(df_plt$Type)) +
      scale_color_brewer(palette = "Reds", direction = -1) +
      theme_my_pub()
  }else{
    return(
      ggplot(df_plt, aes(x = x, y = y, group = Type)) +
        geom_ribbon(aes(ymin = ci_l, ymax = ci_r), fill = "grey80", color="grey80") +
        geom_abline(slope = 1, intercept = 0, color = "black", size = 1) +
        geom_point(aes(color = Type), size = 0.5) +
        labs(x = bquote(Expected -log[10]~italic((P))),
             y = bquote(Observed -log[10]~italic((P))),
             color = group_title) +
        scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"),
                                      RColorBrewer::brewer.pal(8, "Set1")) ) +
        theme_my_pub()
    )
  }
}
