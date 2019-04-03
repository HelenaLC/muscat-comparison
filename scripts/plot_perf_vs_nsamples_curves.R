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
fns <- snakemake@input$sim
fns <- list.files("data/sim_data/kang", full.names = TRUE)
fns <- grep("de10_ns;", fns, value = TRUE)
gis <- .get_gis(fns)

# tidy results
#snakemake@input$res
fns <- list.files("results/kang/de10_ns", full.names = TRUE)
#fns <- grep("runrep=1", fns, value=TRUE)
df <- .tidy_res(fns, gis) %>% 
    data.table %>% split(
        sorted = TRUE, flatten = FALSE,
        by = c("run_rep", "method", "n_samples"))

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

re <- res[[1]]
perf <- .calc_perf(re, facet = "method")
n_samples <- sort(as.numeric(unique(perf@fdrtpr$method)))
cols <- setNames(brewer.pal(length(n_samples), "RdYlGn"), n_samples)

p <- perf %>% reorder_levels(n_samples) %>% 
    prepare_data_for_plot(colorscheme = cols, incloverall = FALSE) %>% 
    .plot_perf_curves + theme(strip.text = element_text(size = 8, face = "bold"))
p$facet$params$ncol <- 5
ggsave(snakemake@output$fig, p,
    width = 15, height = 7.5, units = "cm", 
    dpi = 300, useDingbats = FALSE)


