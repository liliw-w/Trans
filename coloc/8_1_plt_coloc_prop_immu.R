
######################################################
########## 2. visualize the proportions ##########
######################################################

# order the traits based on the number of corresponding regions
#res_coloc_reg_prop = rename(res_coloc_reg_prop, c("Phenocode" = "trait", "trait" = "Phenocode"))

n_region_pval_all_trait = 255

res_coloc_reg_prop = fread(file_res_coloc_reg_prop)
res_coloc_reg_prop$nRegionPval_all_trait = n_region_pval_all_trait

# figure 1: draw bar plot on number of reigons
dat_fig_bar_prop = res_coloc_reg_prop %>%
  pivot_longer(c(nRegion, nRegionPval), names_to = "regionType", values_to = "n") %>%
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc") %>%
  arrange(desc(propPvalColoc), desc(propColoc))
dat_fig_bar_prop$Phenocode = with(dat_fig_bar_prop,
                                  factor(Phenocode,
                                         levels = unique(Phenocode),
                                         labels = sapply(as.character(unique(Phenocode)),
                                                         function(x) strsplit(x, "_", fixed = TRUE)[[1]][1] ))
)

### Don't draw blue regions, only draw green regions, i.e. regions whose lead-SNPs are trans-eQTLs
dat_fig_bar_prop <- dat_fig_bar_prop %>% filter(regionType == "nRegionPval" & regionTypeColoc == "nRegionPvalColoc")


### Draw bar plot (with dual y-axis)
y_lim <- n_region_pval_all_trait
fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = Phenocode)) +
  geom_bar(aes(y = n_region_pval_all_trait), stat = "identity", fill = "#edf5f9") +
  geom_bar(aes(y = n, fill = regionType), stat = "identity", position=position_dodge()) +
  geom_bar(aes(y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = NULL, y = "Regions") +
  scale_fill_brewer(palette = "Paired") + #c("#A6CEE3", "#1F78B4"), c("#B2DF8A", "#33A02C")
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
        axis.text.x = element_text(angle = 60, hjust=1, vjust = 1, size = 10),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.y = element_text(angle=90,vjust =2, size=14),
        axis.title.x = element_text(vjust = -0.2, size=12),
        axis.title.y.right = element_text(angle = 90) )


# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = res_coloc_reg_prop %>%
  select(c(Phenocode, trait, propColoc, propPvalColoc)) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")
dat_fig_line_prop$Phenocode = with(dat_fig_line_prop,
                                  factor(Phenocode,
                                         levels = unique(Phenocode),
                                         labels = sapply(as.character(unique(Phenocode)),
                                                         function(x) strsplit(x, "_", fixed = TRUE)[[1]][1] ))
)


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
saveRDS(fig_combined, 'fig4b_immun.rds')

ggsave(file_plot, fig_combined, width = 6, height = 5)
