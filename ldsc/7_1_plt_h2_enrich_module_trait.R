library(ComplexHeatmap)

file_coloc_fig_order <- '/project2/xuanyao/llw/ldsc/plots/coloc_m_trait_order.rds'
coloc_fig_order <- readRDS(file_coloc_fig_order)


h2_mat <- h2_enrich %>%
  mutate(`Trait Abbreviation` = paste(`Trait Abbreviation`, trait_id, sep = "-")) %>%
  filter(module %in% as.numeric(levels(coloc_fig_order$Module))) %>%
  pivot_wider(
    id_cols = module,
    names_from = `Trait Abbreviation`,
    values_from = Enrichment_p
  ) %>%
  column_to_rownames(var = "module") %>%
  as.matrix() %>%
  t()



col_break <- c(0, 1:4, 10)

annot_trait <- h2_enrich %>%
  mutate(`Trait Abbreviation` = paste(`Trait Abbreviation`, trait_id, sep = "-")) %>%
  distinct(`GWAS Group`, `Trait Abbreviation`) %>%
  left_join(
    distinct(coloc_fig_order, trait_type, trait_color),
    by = c("GWAS Group" = "trait_type")
  ) %>%
  column_to_rownames(var = "Trait Abbreviation")
annot_trait <- annot_trait[rownames(h2_mat), ]


row_ha <- HeatmapAnnotation('Group' = annot_trait$`GWAS Group`,
                            col = list("Group" = setNames(annot_trait$trait_color, annot_trait$`GWAS Group`)),
                            which = "row")


pdf("h2_enrich_order_by_dist.pdf", width = 8.5, height = 5)
Heatmap(-log10(h2_mat), 
        name = "-Log10p", #title of legend
        column_title = "Module", row_title = "Trait",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        column_names_gp = gpar(fontsize = 7),
        right_annotation = row_ha,
        col = circlize::colorRamp2(col_break,
                                   c("white", RColorBrewer::brewer.pal(n = length(col_break), name = "Blues")[-1]))
)
dev.off()

