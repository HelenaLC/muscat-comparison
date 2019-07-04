source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

groups <- c("E <= 0.5", "0.5 < E <= 1", "1 < E <= 2", "E > 2")
.get_group <- function(u) sapply(u, function(v) 
    if (v <= 0.5) 1 else if (v <= 1) 2 else if (v <= 2) 3 else 4) %>% 
    factor

#fns <- list.files("/Users/helena/Dropbox/portmac/results/kang", "ds10;", full.names = TRUE)
#rds <- lapply(fns, readRDS)
res <- lapply(snakemake@input$res, readRDS) %>% 
    map("tbl") %>%
    map(mutate_if, is.factor, as.character) %>% 
    bind_rows %>% 
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
    dplyr::filter(splitval != "overall") %>% 
    mutate_at("splitval", function(u) gsub("group:", "", u)) %>% 
    mutate_at("splitval", factor, labels = groups) %>% 
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) 

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)
# ggsave("/Users/helena/Dropbox/portmac/perf_by_expr.pdf", p,
#     width = 15, height = 6.25, units = "cm", dpi = 300)
saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 6.25,
    dpi = 300, useDingbats = FALSE)