###############################################################
########### Bar plot enrichment in module genes ###########
###############################################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggnewscale)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
fdr_level <- 0.05
n_term_show <- 30
neg_log10_p_show <- 15

module_seq <- list.files(
  '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/enrich_single_module/',
  pattern = 'module_enrich_M\\d+.txt'
) |>
  str_extract_all('\\d+') |> as.numeric()


for (module in module_seq) {
  file_enrich <- paste0('/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/enrich_single_module/module_enrich_M', module, '.txt')
  
  ## output -----
  file_plt_enrich <- paste0('/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/plot/M', module, '_bar.pdf')
  
  
  # read files -----
  enrich <- fread(file_enrich)
  thre_sig <- -log10(fdr_level)
  
  
  # Prepare plot: reorder rows and axis labels -----
  enrich <- mutate(enrich, negative_log10_of_adjusted_p_value = -log10(adjusted_p_value)) %>%
    group_by(gene_module, source) %>%
    slice_max(negative_log10_of_adjusted_p_value, n = n_term_show) %>%
    arrange(negative_log10_of_adjusted_p_value, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(
      term_name = make.unique(term_name), 
      negative_log10_of_adjusted_p_value = sapply(negative_log10_of_adjusted_p_value, min, neg_log10_p_show)
    )
  
  enrich$term_name <- factor(enrich$term_name, levels = enrich$term_name, labels = enrich$term_name)
  enrich <- mutate(enrich, term_name = str_pad(term_name, 100, pad = "."))
  
  
  # plot enrichment of terms -----
  enrich %>%
    ggplot(aes(y = term_name, x = negative_log10_of_adjusted_p_value, fill = negative_log10_of_adjusted_p_value)) +
    geom_col(
      data = subset(enrich, source == "GO:BP"),
      aes(y = term_name, x = negative_log10_of_adjusted_p_value, fill = negative_log10_of_adjusted_p_value),
      show.legend = FALSE, position = position_dodge()
    ) +
    scale_fill_gradient(low = "#c77f7f", high = "#6f0000") +
    
    new_scale_fill() +
    
    geom_col(
      data = subset(enrich, source == "GO:MF"),
      aes(y = term_name, x = negative_log10_of_adjusted_p_value, fill = negative_log10_of_adjusted_p_value),
      show.legend = FALSE, position = position_dodge()
    ) +
    scale_fill_gradient(low = "#d09d7a", high = "#a46934") +
    
    new_scale_fill() +
    
    geom_col(
      data = subset(enrich, source == "KEGG"),
      aes(y = term_name, x = negative_log10_of_adjusted_p_value, fill = negative_log10_of_adjusted_p_value),
      show.legend = FALSE, position = position_dodge()
    ) +
    scale_fill_gradient(low = "#7f99b2", high = "#003366") +
    
    new_scale_fill() +
    
    geom_col(
      data = subset(enrich, source == "REAC"),
      aes(y = term_name, x = negative_log10_of_adjusted_p_value, fill = negative_log10_of_adjusted_p_value),
      show.legend = FALSE, position = position_dodge()
    ) +
    scale_fill_gradient(low = "#bf95cc", high = "#800080") +
    
    geom_vline(xintercept = thre_sig, linetype = "dashed", color = "black") +
    labs(
      y = NULL, 
      x = quote(-Log[10](P[Adjust])), 
    ) +
    scale_y_discrete(
      label = function(x) str_wrap(x, width = 100) |> str_to_title(), limits = enrich$term_name
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_my_pub(legend.position = "right") +
    theme(
      axis.ticks = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 12)
    )
  
  
  # print out key message or write out -----
  ggsave(
    file_plt_enrich,
    width = 5 + 0.25*max(enrich$negative_log10_of_adjusted_p_value), 
    height = max(0.14*nrow(enrich), 2)
  )
  
}

