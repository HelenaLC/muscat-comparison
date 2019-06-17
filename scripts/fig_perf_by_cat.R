source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(purrr)
})

p1 <- readRDS(snakemake@input$gg1)
p2 <- readRDS(snakemake@input$gg2)

df <- map(list(p1, p2), "data") %>% bind_rows(.id = "p_adj") %>% 
    mutate_at("p_adj", factor, labels = paste0("p_adj.", c("loc", "glb")))

p <- .plot_perf_points(df, facet = c("p_adj", "splitval"))
p$facet$params$ncol <- 4

ggsave(snakemake@output$fig, p,
    width = 15, height = 12, units = "cm",
    dpi = 300, useDingbats = FALSE)
