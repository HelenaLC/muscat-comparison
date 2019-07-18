source(snakemake@config$utils)
x <- snakemake@wildcards$x

suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

#fns <- list.files("results/kang", "ds10_ns;", full.names = TRUE)
res <- .read_res(snakemake@input$res) %>% 
    dplyr::mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% setDT %>% 
    split(by = "j", flatten = FALSE) %>% 
    map(group_by, mid) %>% map(function(u) 
        setNames(group_split(u), group_keys(u)[[1]]))

cd <- lapply(seq_along(res), function(i) {
    truth <- res[[i]][[1]][, c(x, "is_de")] %>% 
        data.frame(row.names = NULL, check.names = FALSE)
    pvals <- lapply(c("p_val", "p_adj.loc"), map, .x = res[[i]]) %>% 
        map(data.frame, check.names = FALSE)
    dfs <- c(list(truth), pvals)
    names(dfs) <- c("truth", "pval", "padj")
    do.call(COBRAData, dfs)
})

perf <- lapply(cd, calculate_performance,
    binary_truth = "is_de", 
    aspects = c("fdrtpr", "fdrtprcurve"),
    splv = x, maxsplit = Inf)

df <- map(perf, "fdrtpr") %>% 
    bind_rows(.id = "j") %>% 
    dplyr::select(splitval, thr, method, TPR, FDR) %>% 
    dplyr::filter(splitval != "overall") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>%  
    mutate_at("splitval", function(u) {
        u <- gsub(paste0(x, ":"), "", u)
        v <- sort(unique(as.numeric(u)))
        factor(u, levels = v)
    }) %>% 
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)

saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    width = 15, height = 6.2, units = "cm",
    dpi = 300, useDingbats = FALSE)
