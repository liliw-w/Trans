ggplot() +
geom_violin() +
geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") +
stat_summary(fun.data = mean_sdl, geom = "pointrange") +
geom_quasirandom(varwidth = TRUE, shape = 16, size = 1, alpha = 0.7)

