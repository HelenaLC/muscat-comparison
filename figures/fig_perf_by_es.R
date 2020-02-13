suppressMessages({
    library(cowplot)
    library(ggplot2)
})

fns <- paste0(config$dids, "-perf_by_es_loc.rds")
ps <- lapply(file.path("plots", fns), readRDS)
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(ncol = 1,
    plotlist = c(ps, list(lgd)),
    rel_heights = c(2, 2, 0.2),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold",
    align = "h", axis = "lr")

ggsave(file.path("figures", "perf_by_es.pdf"), 
    p, dpi = 300, useDingbats = FALSE,
    width = 15, height = 30, units = "cm")
