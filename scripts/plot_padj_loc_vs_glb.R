# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(ggplot2)
        library(purrr)
    }))

# source utils
source(snakemake@input$utils)

# load truth & tidy results
gis <- .get_gis(snakemake@input$sim)
res <- .tidy_res(snakemake@input$res, gis) %>% 
    setDT %>% split(
        sorted = TRUE, flatten = FALSE,
        by = c("sim_rep", "method"))

# caluclate performance for each
# replicate, method, and p-adjustment type
ps <- paste0("p_adj.", c("loc", "glb"))
names(ps) <- ps
perf <- lapply(ps, function(p) {
    ps <- c("p_val", p)
    names(ps) <- ps
    lapply(seq_along(res), function(i)
        lapply(ps, function(p)
            res[[i]] %>% 
                map_depth(1, pull, p) %>% 
                bind_cols) %>% 
            c(list(truth = data.frame(
                is_de = gis[[i]]$is_de,
                rep = i)))) %>% 
        map_depth(2, data.frame, 
            check.names = FALSE,
            stringsAsFactors = FALSE) %>% 
        map(.calc_perf)
})

# prep. data.frame for plotting
df <- map_depth(perf, 2, function(u)
    u %>% fdrtpr %>%
        select(method, thr, FDR, TPR) %>% 
        mutate_if(is.factor, as.character)) %>% 
    map_depth(1, bind_rows, .id = "rep") %>% 
    bind_rows(.id = "p_adj") %>% 
    mutate_at("thr", function(u)
        as.numeric(gsub("thr", "", u))) %>% 
    # average over replicates
    group_by(method, thr, p_adj) %>% 
    summarize_at(c("FDR", "TPR"), mean) %>% ungroup %>% 
    mutate_at("method", factor, levels = names(method_colors))

p <- .plot_perf_points(df, facet = "p_adj")
p$facet$params$ncol <- 2

ggsave(snakemake@output$fig, p,
    width = 10, height = 8, units = "cm",
    dpi = 300, useDingbats = FALSE)
