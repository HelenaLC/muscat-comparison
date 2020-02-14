suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

# wcs = list(x = "c", did = "magl")
# args <- list(
#     res = list.files("results", sprintf("%s,de10_n%s,", wcs$did, wcs$x), full.names = TRUE),
#     ggp = file.path("plots", paste0(wcs$did, sprintf("-perf_by_n%s.rds", wcs$x))),
#     fig = file.path("plots", paste0(wcs$did, sprintf("-perf_by_n%s.pdf", wcs$x))))

res <- .read_res(args$res) %>% 
    dplyr::mutate(id = sprintf(
        "i%sj%sc%ss%s%s%s", 
        i, j, c, s, gene, cluster_id)) %>% 
    dplyr::mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% setDT %>% 
    split(by = "i", flatten = FALSE) %>% 
    map(group_by, mid) %>% map(function(u) 
        setNames(group_split(u), group_keys(u)[[1]]))

# some methods may fail for too low number of cells / replicates;
# the below chunk fills in missing results to match in dimension
for (i in seq_along(res)) {
    u <- res[[i]][[1]]
    any_missing <- which(sapply(res[[i]], nrow) != nrow(u))
    for (j in any_missing) {
        v <- res[[i]][[j]]
        m <- match(setdiff(u$id, v$id), u$id)
        filler <- mutate_if(u[m, ], is.numeric, 
            function(u) replace(u, TRUE, NA))
        v <- rbind(v, filler)
        v <- v[match(u$id, v$id), ]
        res[[i]][[j]] <- v
    }
}

cd <- lapply(seq_along(res), function(i) {
    truth <- res[[i]][[1]][, c("id", "is_de", wcs$x)] %>% 
        data.frame(row.names = NULL, check.names = FALSE)
    ps <- lapply(c("p_val", "p_adj.loc"), map, .x = res[[i]]) %>% 
        map(data.frame, check.names = FALSE)
    dfs <- c(list(truth), ps)
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
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)
