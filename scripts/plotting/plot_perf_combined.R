# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(cowplot)
        library(ggplot2)
    }))

# load ggplot-objects
fns <- list.files("figures/kang", ".rds", full.names = TRUE)[c(1, 3)]
ps <- lapply(fns, readRDS) 

# arrange plots
lgd <- get_legend(ps[[1]] + theme(
    legend.spacing = unit(rep(0, 4), "mm")))
ps <- lapply(ps, function(p) 
    p + theme(legend.position = "none",
        plot.margin = unit(c(0,0,2,4), "mm")))
ps <- c(ps, list(lgd))

p <- plot_grid(plotlist = ps,
    ncol = 1, axis = "t", align = "v",
    rel_heights = c(1.8, 1, 0.4),
    labels = c("a", "b", ""),
    label_size = 10)

# save figure to .pdf
ggsave("figures/kang/perf_combined.pdf", p,
    width = 15, height = 16, units = "cm",
    dpi = 300, useDingbats = FALSE)

# ps <- lapply(fns, readRDS) 
# lgd <- get_legend(
#     ps[[1]] + guides(color = guide_legend(ncol = 1)) +
#         theme(legend.spacing = unit(rep(0, 4), "mm")))
# ps <- lapply(ps, function(p) 
#     p + theme(legend.position = "none",
#         plot.margin = unit(c(0,0,2,4), "mm")))
# ps[[2]]$facet$params$ncol <- 2
# ps[[3]]$facet$params$ncol <- 2
# 
# p1 <- plot_grid(plotlist = c(ps[1], list(lgd)), 
#     rel_widths = c(3, 1), axis = "tl", align = "h",
#     labels = "a",
#     label_size = 10)
# p2 <- plot_grid(plotlist = ps[2:3], 
#     labels = c("a", "b"),
#     label_size = 10)
# p <- plot_grid(p1, p2, ncol = 1,
#     axis = "tblr", align = "hv",
#     rel_heights = c(1, 1))
# ggsave("figures/kang/perf_combine2.pdf", p,
#     width = 16, height = 14, units = "cm",
#     dpi = 300, useDingbats = FALSE)
