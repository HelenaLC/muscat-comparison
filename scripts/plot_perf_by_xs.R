suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(magrittr)
    library(purrr)
})

#args <- list(res = list.files("~/projects/portmac/results", "kang,de10_gs[0-9],", full.names = TRUE))
res <- .read_res(args$res, wcs$inc) %>% 
    dplyr::mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% setDT %>% 
    split(by = c("sid", "i", "mid"), flatten = FALSE)

cd <- map_depth(res, 2, function(u) {
    truth <- setDF(select(u[[1]], c("sid", "i", "is_de")))
    pvals <- lapply(c("p_val", "p_adj.loc"), map, .x = u)
    rmv <- map_depth(pvals, 1, vapply, is.null, logical(1))
    pvals <- map_depth(pvals, 1, function(u) 
        setDF(u[!vapply(u, is.null, logical(1))]))
    dfs <- set_names(c(list(truth), pvals), c("truth", "pval", "padj"))
    do.call(COBRAData, dfs)
})

perf <- map_depth(cd, 2, calculate_performance,
    aspects = "fdrtpr", binary_truth = "is_de")

gg_df <- map_depth(perf, 2, fdrtpr) %>% 
    map(bind_rows, .id = "i") %>% 
    bind_rows(.id = "sid") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>% 
    group_by(sid, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

p <- .plot_perf_points(gg_df, facet = "sid")
p$facet$params$ncol <- nlevels(factor(gg_df$sid))

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)