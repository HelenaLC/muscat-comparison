suppressMessages({
    library(cowplot)
    library(ggplot2)
})

ps <- lapply(args$ggp, readRDS)
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(
    plotlist = c(ps, list(lgd)),
    ncol = 1, align = "v", axis = "t",
    rel_heights = c(5, 5, 1),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

ggsave(args$fig, p,
    width = 15, height = 14, units = "cm",
    dpi = 300, useDingbats = FALSE)
