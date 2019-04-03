# plot method runtimes vs. nb. of genes
# ==============================================================================

# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(purrr)
    })
)

# source utils
source(snakemake@input$utils)

rds <- list.files(snakemake@input$res, full.names = TRUE)
fns <- basename(rds)
res <- lapply(rds, readRDS)
res <- map(res, "rt")
res <- map(res, sum)
rmv <- vapply(res, is.null, logical(1))
fns <- fns[!rmv]
res <- res[!rmv]

zeallot::`%<-%`(
    c(method, n_genes, n_cells, rep), 
    .parse_fns(fns))

df <- data.frame(
    time = res %>% unlist,
    method, n_genes, n_cells, rep)

df$method <- factor(df$method, levels = names(method_colors))
df$id <- with(df, paste0(method, n_genes))

# plotting ---------------------------------------------------------------------

p <- ggplot(df, aes(x = n_genes, y = time, group = id, col = method, fill = method)) +
    geom_boxplot(size = 0.2, outlier.size = 0.2, alpha = 0.4) + 
    stat_smooth(formula = y ~ x, alpha = 0.1, size = 0.2, aes(group = method)) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    scale_x_continuous(breaks = unique(n_genes), trans = "log2", expand = c(0,.04)) + 
    scale_y_continuous(breaks = c(1,10,30,60,120,240,480), trans = "log10", expand = c(0,.06)) +
    labs(x = "nb. of genes", y = "time (s)") +
    prettify() + theme(aspect.ratio = 2/3)

ggsave(snakemake@output$fig, p,
    width = 16, height = 8, units = "cm", 
    dpi = 300, useDingbats = FALSE)    
