suppressMessages({
    library(cowplot)
    library(ggplot2)
})

ref_id <- "kang"

ps <- list(
    readRDS(file.path("figures", paste0(ref_id, "-perf_by_cat.rds"))),
    readRDS(file.path("plots", paste0(ref_id, "-perf_by_nc.rds"))))

lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(ncol = 1, 
    plotlist = c(ps, list(lgd)),
    rel_heights = c(5.4, 3, 0.9),
    labels = c("a", "b", ""),
    label_size = 10,
    label_fontface = "bold")

fn <- file.path("figures", paste0(ref_id, "-perf_combined.pdf"))
ggsave(fn, p,
    width = 15, height = 15, units = "cm", 
    dpi = 300, useDingbats = FALSE)   






