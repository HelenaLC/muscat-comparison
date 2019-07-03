source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

#lfc_groups <- c("|logFC| < 0.5", "0.5 <= |logFC| < 1", "", "1 >= |logFC| < 2", "|logFC| >= 2")
.get_group <- function(u) sapply(u, function(v) 
    if (v < 0.2) 1 else if (v < 0.4) 2 else if (v < 0.6) 3 else if (v < 0.8) 4 else if (v < 1) 5 else 6) %>% 
    factor(levels = seq_len(7) - 1)

#fns <- list.files("/Users/helena/Dropbox/portmac/results/magl", "ds10;", full.names = TRUE)
res0 <- rds %>% 
    map("tbl") %>% bind_rows %>% 
    mutate(group = .get_group(abs(.$sim_lfc))) %>% 
    {.$group[!.$is_de] <- 0; return(.)} %>% 
    setDT %>% split(by = c("i", "group", "mid"), flatten = FALSE)

sub <- map_depth(res0, 3, function(u) {
    i <- as.character(u$i[1])
    m <- as.character(u$mid[1])
    v <- bind_rows(u, res0[[i]][["0"]][[m]])
    v$group <- u$group[1]
    return(v)
})
sub <- map(sub, `[`, -1)
res <- sub
#res <- map(res0, `[`, -1)
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
    dplyr::filter(splitval != "overall") %>%
    group_by(splitval, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) 

p <- .plot_perf_points(df)
p$facet$params$ncol <- nlevels(df$splitval)
ggsave("/Users/helena/Dropbox/portmac/perf_by_lfc.pdf", p,
    width = 15, height = 12, units = "cm", dpi = 300)

saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)
