# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(ggplot2)
        library(purrr)
        library(reshape2)
        library(tibble)
    })
)

# source utils
source(snakemake@input$utils)

# load truth & tidy results
gis <- .get_gis(snakemake@input$sim)
df <- .tidy_res(snakemake@input$res, gis) %>% 
    data.table %>% split(
        sorted = TRUE, flatten = FALSE,
        by = c("run_rep", "method", "n_cells"))

# prep. data for iCOBRA
ps <- c("p_val", "p_adj.loc")
names(ps) <- c("p_val", "p_adj")
truth_cols <- c("is_de", "method")
res <- lapply(seq_along(df), function(i) {
    re <- lapply(ps, function(p) {
        df[[i]] %>% map_depth(2, pull, p) %>% 
            map_depth(1, bind_cols) %>% 
            map_depth(1, add_column, is_de = gis[[1]]$is_de) %>% 
            bind_rows(.id = "method")
    })
    list(truth = select(re[[1]], truth_cols)) %>% 
        c(map(re, select, -truth_cols))
}) %>% map_depth(2, data.frame, check.names = FALSE)

# calculate performances
df <- lapply(seq_along(res), function(i) {
    perf <- .calc_perf(res[[i]], facet = "method", thrs = 0.05)
    fdrtpr(perf) %>% select(c("TPR", "FDR", "method", "splitval")) %>% 
        rename(method = "n_cells", splitval = "method") %>% 
        filter(method != "overall") %>% 
        mutate(rep = i)
}) %>% bind_rows

# prep. for plotting
gg_df <- df %>% 
    mutate_at("n_cells", as.numeric) %>% 
    mutate_at("method", function(u) 
        gsub("^method:", "", u) %>% 
            factor(levels = names(method_colors))) %>% 
    melt(id.vars = c("method", "n_cells", "rep"))

p <- ggplot(gg_df, aes(x = n_cells, y = value, 
    group = paste(method, n_cells), col = method, fill = method)) + 
    facet_wrap(~ variable, scales = "free", ncol = 3) + 
    geom_boxplot(width = 0.8, alpha = 0.8, size = 0.2, outlier.size = 0.1) + 
    scale_color_manual(NULL, values = method_colors) +
    scale_fill_manual(NULL, values = method_colors) + 
    guides(col = guide_legend(ncol = 3)) +
    scale_x_continuous(trans = "log2", breaks = unique(gg_df$n_cells)) +
    scale_y_continuous(NULL, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    prettify() + theme(aspect.ratio = 2/3,
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(snakemake@output$fig, p,
    width = 15, height = 8, units = "cm", 
    dpi = 300, useDingbats = FALSE)
