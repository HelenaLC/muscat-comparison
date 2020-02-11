suppressMessages({
    library(cowplot)
    library(ggplot2)
})

fns <- paste0(config$did, "-sim_vs_est_lfc.rds")
ps <- lapply(file.path("plots", fns), readRDS)
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(
    legend.position = "none", 
    strip.text = element_text(size = 4)))

ps[[1]]$layers[[1]]$geom_params$raster.dpi <- 50
ps[[2]]$layers[[1]]$geom_params$raster.dpi <- 50

p <- plot_grid(ncol = 1,
    plotlist = c(ps, list(lgd)),
    rel_heights = c(15, 15, 1),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

ggsave(file.path("figures", "sim_vs_est_lfc.pdf"), 
    p, width = 15, height = 20, units = "cm",
    dpi = 200, useDingbats = FALSE)
