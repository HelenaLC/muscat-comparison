suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

#args <- list(res = list.files("~/projects/portmac/results", "kang,de10_ns", full.names = TRUE))
res <- .read_res(args$res) %>% 
    dplyr::mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% setDT %>% 
    split(by = "j", flatten = FALSE) %>% 
    map(group_by, mid) %>% map(function(u) 
        setNames(group_split(u), group_keys(u)[[1]]))

cd <- lapply(seq_along(res), function(i) {
    truth <- res[[i]][[1]][, c(wcs$x, "is_de")] %>% 
        data.frame(row.names = NULL, check.names = FALSE)
    pvals <- lapply(c("p_val", "p_adj.loc"), map, .x = res[[i]]) %>% 
        map(data.frame, check.names = FALSE)
    dfs <- c(list(truth), pvals)
    names(dfs) <- c("truth", "pval", "padj")
    do.call(COBRAData, dfs)
})

perf <- lapply(cd, calculate_performance, 
    aspects = "fdrtpr", binary_truth = "is_de", 
    splv = wcs$x, maxsplit = Inf)

df <- map(perf, "fdrtpr") %>% 
    bind_rows(.id = "j") %>% 
    dplyr::select(splitval, thr, method, TPR, FDR) %>% 
    dplyr::filter(splitval != "overall") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>%  
    mutate_at("splitval", function(u) {
        u <- gsub(paste0(wcs$x, ":"), "", u)
        v <- sort(unique(as.numeric(u)))
        factor(u, levels = v)
    }) %>% 
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 6.2, units = "cm",
    dpi = 300, useDingbats = FALSE)
