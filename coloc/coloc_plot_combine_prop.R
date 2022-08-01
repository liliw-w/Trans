######################################################
########## 2. visualize the proportions ##########
######################################################
pvalThre <- 'module_QTL_sig'
file_plot <- 'ukbb_immun_prop.pdf'
n_region_pval_all_trait <- 255


file_res_coloc_reg_prop_b <- paste0("/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/coloc_region_prop_pvalThre-", pvalThre, ".txt")
file_res_coloc_reg_prop_immun <- paste0("~/xuanyao_llw/coloc/immune_traits/pmid_all/coloc_region_prop_pvalThre-", pvalThre, ".txt")


res_coloc_reg_prop_b <- fread(file_res_coloc_reg_prop_b)
res_coloc_reg_prop_immun <- fread(file_res_coloc_reg_prop_immun)


res_coloc_reg_prop_b <- mutate(res_coloc_reg_prop_b,
                               "Phenocode" = trait,
                               "trait_type" = "immune")
res_coloc_reg_prop_immun <- mutate(res_coloc_reg_prop_immun,
                                   "Phenocode" = sapply(Phenocode,
                                                        function(x) strsplit(x, "_", fixed = TRUE)[[1]][1]),
                                   "trait_type" = "blood"
                                   )


# order the traits based on the number of corresponding regions
res_coloc_reg_prop <- rbind(res_coloc_reg_prop_b, res_coloc_reg_prop_immun)
res_coloc_reg_prop$trait_type <- factor(res_coloc_reg_prop$trait_type,
                                        levels = c("immune", "blood"))
res_coloc_reg_prop$nRegionPval_all_trait <- n_region_pval_all_trait



dat_fig_bar_prop <- res_coloc_reg_prop %>%
  pivot_longer(c(nRegion, nRegionPval), names_to = "regionType", values_to = "n") %>%
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc") %>%
  group_by(trait_type) %>%
  arrange(desc(propPvalColoc), desc(propColoc), .by_group = TRUE) %>%
  ungroup()
dat_fig_bar_prop$Phenocode <- factor(dat_fig_bar_prop$Phenocode,
                                     levels = unique(dat_fig_bar_prop$Phenocode))


### Don't draw blue regions, only draw green regions, i.e. regions whose lead-SNPs are trans-eQTLs
dat_fig_bar_prop <- dat_fig_bar_prop %>% filter(regionType == "nRegionPval" & regionTypeColoc == "nRegionPvalColoc")


### Draw bar plot (with dual y-axis)
y_lim <- n_region_pval_all_trait
fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = Phenocode)) +
  geom_bar(aes(y = n_region_pval_all_trait), stat = "identity", fill = "#edf5f9") + #"#eaf4e5"
  geom_bar(aes(y = n, fill = regionType), stat = "identity", position=position_dodge()) +
  geom_bar(aes(y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = NULL, y = "Regions") +
  scale_fill_brewer(palette = "Paired") +
  #scale_fill_manual(values = c("#99cc7f", "#2d8900")) +
  scale_y_continuous(limits = c(0, y_lim),
                     sec.axis = sec_axis(~./y_lim, name = "Coloc Proportion")
  ) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour="black"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text.x = element_text(angle = 60, hjust=1, vjust = 1, size = 8),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(angle = 90,vjust = 2, size = 12),
        axis.title.x = element_text(vjust = -0.2, size = 12),
        axis.title.y.right = element_text(angle = 90) )


# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop <- res_coloc_reg_prop %>%
  select(c(Phenocode, trait, trait_type, propColoc, propPvalColoc)) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")


### Don't draw blue regions, only draw green regions, i.e. regions whose lead-SNPs are trans-eQTLs
dat_fig_line_prop <- dat_fig_line_prop %>% filter(Type == "propPvalColoc")


### Add the proportion line on the second y-axis of the above bar plot
fig_combined <- fig_bar_prop +
  geom_line(data = dat_fig_line_prop,
            aes(x = Phenocode,
                y = proportion*y_lim, group = Type), color = "blue") +
  geom_point(data = dat_fig_line_prop,
             aes(x = Phenocode,
                 y = proportion*y_lim, group = Type), color = "blue") +
  theme(axis.text.x = element_text(color = "black"))


### save figure object and figure
saveRDS(fig_combined, 'ukbb_immun_prop.rds')

ggsave(file_plot, fig_combined, width = 7, height = 3)
