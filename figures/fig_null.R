suppressMessages({
    library(cowplot)
    library(ggplot2)
})

fns <- paste0(config$dids, "-null.rds")
ps <- lapply(file.path("plots", fns), readRDS)
lgd <- get_legend(ps[[1]] + guides(col = FALSE, 
    lty = guide_legend(ncol = 1, order = 2),
    fill = guide_legend(ncol = 4, order = 1,
        override.aes = list(alpha = 1, col = NA))) +
        theme(legend.position = "bottom"))
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(
    plotlist = c(ps, list(lgd)),
    ncol = 1, align = "v", axis = "t",
    rel_heights = c(5, 5, 0.8),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

ggsave(file.path("figures", "null.pdf"), p,
    width = 15, height = 20, units = "cm",
    dpi = 300, useDingbats = FALSE)
