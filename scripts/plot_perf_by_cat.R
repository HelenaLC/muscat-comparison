suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

# wcs <- list(padj = "loc", inc = "treat")
# args <- list(res = list.files("results", "kang,d[a-z][0-9]+,", full.names = TRUE))
res <- .read_res(args$res, wcs$inc) %>%  
    mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% setDT %>% 
    split(by = c("i", "sid", "mid"), flatten = FALSE)

p_adj <- paste0("p_adj.", wcs$padj)
cd <- lapply(seq_along(res), function(i) COBRAData( 
    pval = as.data.frame(bind_rows(map(map_depth(res[[i]], 2, "p_val"), bind_cols))),
    padj = as.data.frame(bind_rows(map(map_depth(res[[i]], 2, p_adj), bind_cols))),
    truth = data.frame(
        row.names = NULL,
        sim_id = unlist(map(map_depth(res[[i]], 2, "sid"), 1)),
        is_de = unlist(map(map_depth(res[[i]], 2, "is_de"), 1)))))

perf <- lapply(cd, calculate_performance, 
    aspects = "fdrtpr", binary_truth = "is_de", 
    splv = "sim_id", maxsplit = Inf)

df <- map(perf, function(u) 
    select(fdrtpr(u), splitval, thr, method, TPR, FDR)) %>% 
    bind_rows(.id = "i") %>% 
    mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>% 
    mutate_at("method", factor, levels = names(.meth_cols)) %>% 
    mutate_at("splitval", function(u) {
        u <- gsub("sim_id:([a-z]+)[0-9]+", "\\1", u)
        factor(u, 
            levels = c("de", "dp", "dm", "db"),
            labels = c("DE", "DP", "DM", "DB"))
    }) %>% 
    dplyr::filter(splitval != "overall") %>%
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean)

p <- .plot_perf_points(df, include = wcs$inc)
p$facet$params$ncol <- nlevels(df$splitval)

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)