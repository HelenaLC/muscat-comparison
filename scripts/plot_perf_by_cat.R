source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

#fns <- list.files("/users/helena/dropbox/portmac/results/kang", "d[a-z][0-9]+;", full.names = TRUE)
res <- lapply(snakemake@input$res, readRDS) %>% 
    map("tbl") %>% 
    map(mutate_if, is.factor, as.character) %>% 
    bind_rows %>% setDT %>% 
    split(by = c("i", "sid", "mid"), flatten = FALSE)

p_adj <- paste0("p_adj.", snakemake@wildcards$padj)

cd <- lapply(seq_along(res), function(i) COBRAData( 
    pval = as.data.frame(bind_rows(map(map_depth(res[[i]], 2, "p_val"), bind_cols))),
    padj = as.data.frame(bind_rows(map(map_depth(res[[i]], 2, p_adj), bind_cols))),
    truth = data.frame(
        row.names = NULL,
        sim_id = unlist(map(map_depth(res[[i]], 2, "sid"), 1)),
        is_de = unlist(map(map_depth(res[[i]], 2, "is_de"), 1)))))

perf <- lapply(cd, calculate_performance,
    aspects = c("fdrtpr", "fdrtprcurve"),
    splv = "sim_id", maxsplit = Inf,
    binary_truth = "is_de")

df <- map(perf, function(u) 
    select(fdrtpr(u), splitval, thr, method, TPR, FDR)) %>% 
    bind_rows(.id = "i") %>% 
    mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>% 
    mutate_at("method", factor, levels = names(.meth_cols)) %>% 
    mutate_at("splitval", function(u) {
        u <- gsub("sim_id:([a-z]+)[0-9]+", "\\1", u)
        factor(toupper(u), levels = c("DS", "DP", "DM", "DB"))
    }) %>% 
    dplyr::filter(splitval != "overall") %>%
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) 

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)

saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)
