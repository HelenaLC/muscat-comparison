suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

groups <- c("E <= 0.1", "0.1 < E <= 0.5", "0.5 < E <= 1", "E > 1")
.get_group <- function(u) sapply(u, function(v) 
    if (v <= 0.1) 1 else if (v <= 0.5) 2 else if (v <= 1) 3 else 4) %>% 
    factor

#fns <- list.files("results/kang", "d[a-z]10;", full.names = TRUE)
res <- .read_res(args$res) %>% 
    dplyr::mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::mutate(group = .get_group(.$E)) %>% setDT %>% 
    split(by = c("sid", "i", "group", "mid"), flatten = FALSE)

p_adj <- paste0("p_adj.", wcs$padj)
cd <- lapply(names(res), function(sid) {
    lapply(seq_along(res[[sid]]), function(i) COBRAData( 
        pval = as.data.frame(bind_rows(map(map_depth(res[[sid]][[i]], 2, "p_val"), bind_cols))),
        padj = as.data.frame(bind_rows(map(map_depth(res[[sid]][[i]], 2, p_adj), bind_cols))),
        truth = data.frame(
            row.names = NULL,
            group = unlist(map(map_depth(res[[sid]][[i]], 2, "group"), 1)),
            is_de = unlist(map(map_depth(res[[sid]][[i]], 2, "is_de"), 1)))))
})

perf <- map_depth(cd, 2, calculate_performance, 
    aspects = "fdrtpr", binary_truth = "is_de",
    splv = "group", maxsplit = Inf)

df <- map_depth(perf, 2, fdrtpr) %>% 
    map(bind_rows, .id = "i") %>% 
    bind_rows(.id = "sid") %>% 
    select(sid, splitval, thr, method, TPR, FDR) %>% 
    mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>%
    mutate_at("method", factor, levels = names(.meth_cols)) %>% 
    dplyr::filter(splitval != "overall") %>% 
    mutate_at("splitval", function(u) gsub("group:", "", u)) %>% 
    mutate_at("splitval", factor, labels = groups) %>% 
    group_by(splitval, sid, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    ungroup %>% mutate_at("sid", factor, 
        labels = gsub("10", "", names(res))) %>% 
    mutate_at("sid", factor, 
        levels = c("de", "dp", "dm", "db"),
        labels = c("DE", "DP", "DM", "DB"))

p <- .plot_perf_points(df) +
    facet_grid(rows = vars(sid), cols = vars(splitval))

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 16.2, units = "cm",
    dpi = 300, useDingbats = FALSE)
