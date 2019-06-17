source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

x <- snakemake@wildcards$x
#fns <- list.files("results/magl", "ds10_nc;", full.names = TRUE)
res <- lapply(snakemake@input$res, readRDS) %>% 
    map("tbl") %>% bind_rows %>% 
    setDT %>% split(by = "j", flatten = FALSE)

cd <- lapply(seq_along(res), function(i) {
    res2 <- group_by(res[[i]], mid) %>% 
        {set_names(group_split(.), group_keys(.)[[1]])}
    dfs <- c(
        lapply(c("p_val", "p_adj.loc"), map, .x = res2),
        list(truth = res2[[1]][, c("is_de", x)])) %>% 
        map(data.frame, check.names = FALSE)
    do.call(COBRAData, dfs)
})

perf <- lapply(cd, calculate_performance,
    binary_truth = "is_de", 
    aspects = c("fdrtpr", "fdrtprcurve"),
    splv = x, maxsplit = Inf)

df <- map(perf, "fdrtpr") %>% 
    bind_rows(.id = "j") %>% 
    select(splitval, thr, method, TPR, FDR) %>% 
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
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)