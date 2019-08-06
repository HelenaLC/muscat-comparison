config <- yaml::read_yaml("config.yaml")
source(config$utils)

suppressMessages({
    library(cowplot)
    library(DeLuciatoR)
    library(ggplot2)
    library(purrr)
})

ps <- sapply(config$dids, function(id) 
    lapply(file.path("figures", id, "null.rds"), readRDS))
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(
    plotlist = c(ps, list(lgd)),
    ncol = 1, align = "v", axis = "t",
    rel_heights = c(5, 5, 1),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

ggsave(file.path("figures", "null.pdf"), p,
    width = 15, height = 14, units = "cm",
    dpi = 300, useDingbats = FALSE)
