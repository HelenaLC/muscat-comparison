source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

groups <- c("E <= 0.1", "0.1 < E <= 0.2", "0.2 < E <= 0.5", "0.5 < E <= 1", "1 < E <= 2", "E > 2")
.get_group <- function(u) sapply(u, function(v) 
    if (v <= 0.1) 1 else if (v <= 0.2) 2 else if (v <= 0.5) 3 else if (v <= 1) 4 else if (v <= 2) 5 else 6) %>% 
    factor

#fns <- list.files("/Users/helena/Dropbox/portmac/results/magl", "ds10;", full.names = TRUE)
#rds <- lapply(fns, readRDS)
res <- rds %>% 
    map("tbl") %>% bind_rows %>% 
    mutate(avg_expr = sim_mean.A + sim_mean.B) %>% 
    mutate(group = .get_group(abs(.$avg_expr))) %>% 
    setDT %>% split(by = c("i", "group", "mid"), flatten = FALSE)

p_adj <- paste0("p_adj.", snakemake@wildcards$padj)

cd <- lapply(seq_along(res), function(i) COBRAData( 
    pval = as.data.frame(bind_rows(map(map_depth(res[[i]], 2, "p_val"), bind_cols))),
    padj = as.data.frame(bind_rows(map(map_depth(res[[i]], 2, p_adj), bind_cols))),
    truth = data.frame(
        row.names = NULL,
        group = unlist(map(map_depth(res[[i]], 2, "group"), 1)),
        is_de = unlist(map(map_depth(res[[i]], 2, "is_de"), 1)))))

perf <- lapply(cd, calculate_performance,
    aspects = c("fdrtpr", "fdrtprcurve"),
    splv = "group", maxsplit = Inf,
    binary_truth = "is_de")

df <- map(perf, function(u) 
    select(fdrtpr(u), splitval, thr, method, TPR, FDR)) %>% 
    bind_rows(.id = "i") %>% 
    mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>% 
    mutate_at("method", factor, levels = names(.meth_cols)) %>% 
    mutate_at("splitval", function(u) factor(gsub("group:", "", u))) %>% 
        #factor(, levels = groups)) %>% 
    dplyr::filter(splitval != "overall") %>%
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) 

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)
ggsave("/Users/helena/Dropbox/portmac/perf_by_expr.pdf", p,
    width = 24, height = 8, units = "cm", dpi = 300)
saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)