source('~/Trans/plot/theme_my_pub_hist.R')
tmp <- readRDS('~/xuanyao_llw/MODULES/MSigDB/result/coexp.module.rds')

str(tmp)

coexp_module_annot <- enframe(
  coexp_module$moduleName,
  name = "annot_module",
  value = "module"
)


tmp <- tmp$moduleLabels %>%
  enframe(name = "gene", value = "module") %>%
  left_join(
    coexp_module_annot,
    by = "module"
  )%>%
  arrange(module)
tmp$annot_module <- factor(tmp$annot_module, levels = unique(tmp$annot_module))


ggplot(tmp, aes(x = annot_module)) +
  geom_bar() +
  labs(x = "Pathways", y = "Pathway size") +
  theme_my_pub() +
  theme_my_pub_hist +
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.ticks.x = element_blank()
  )


ggsave(
  "/project2/xuanyao/llw/MODULES/MSigDB/plot/module_size.pdf",
  width = 7, height = 6
)

