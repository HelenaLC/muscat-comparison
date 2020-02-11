suppressMessages({
    library(cowplot)
    library(ggplot2)
})

fns <- paste0(config$dids, "-pb_mean_disp.rds")
ps <- lapply(file.path("plots", fns), readRDS)
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(
    legend.position = "none", 
    aspect.ratio = NULL))

ps[[1]]$layers[[1]]$geom_params$raster.dpi <- 50
ps[[2]]$layers[[1]]$geom_params$raster.dpi <- 50

ps[[1]]$layers[[2]]$geom_params$raster.dpi <- 50
ps[[2]]$layers[[2]]$geom_params$raster.dpi <- 50
    
p <- plot_grid(ncol = 1,
    plotlist = c(ps, list(lgd)),
    rel_heights = c(1.5, 1.5, 0.2),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

ggsave(file.path("figures", "pb_mean_disp.pdf"), 
    p, width = 15, height = 8, units = "cm",
    dpi = 300, useDingbats = FALSE)
