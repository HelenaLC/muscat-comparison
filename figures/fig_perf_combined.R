config <- yaml::read_yaml("config.yaml")
source(config$utils)

suppressMessages({
    library(cowplot)
    library(ggplot2)
})

id <- "kang"

p1 <- readRDS(file.path("figures", id, "perf_by_cat.rds"))
p2 <- readRDS(file.path("figures", id, "perf_by_nc.rds"))

ps <- list(p1, p2)
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(ncol = 1, 
    plotlist = c(ps, list(lgd)),
    rel_heights = c(5.4, 3, 0.9),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

ggsave(file.path("figures", id, "perf_combined.pdf"), p,
    width = 15, height = 15, units = "cm", 
    dpi = 300, useDingbats = FALSE)   






