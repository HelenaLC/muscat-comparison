source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(magrittr)
    library(purrr)
})

fns <- list.files(snakemake@config$sim_pars, "ds10_ss[0-9]+;", full.names = TRUE)
sim_pars <- lapply(fns, yaml::read_yaml)
print(sim_pars)
n <- sapply(c("nk", "ns", "nc"), function(u) unlist(map(sim_pars, u)))
n <- with(as.data.frame(n), nc / (2 * nk * ns))
ss <- map(map(sim_pars, "probs"), 2) 
ss <- sapply(seq_along(ss), function(i) n[i] * ss[[i]]) 
colnames(ss) <- gsub(".json", "", basename(fns))

gg_df <- melt(ss) %>% 
    set_colnames(c("sample_id", "sim_id", "n")) %>% 
    replicate(n = 2, simplify = FALSE) %>% 
    bind_rows(.id = "group_id") %>% 
    mutate_at("group_id", function(u) factor(u, labels = c("A", "B"))) %>% 
    mutate(sample_id = factor(paste0(group_id, sample_id)))

hists <- ggplot(gg_df, aes(x = sample_id, y = n, fill = sample_id)) +
    facet_wrap("sim_id", nrow = 1) +
    geom_bar(size = 1, width = 0.5, stat = "identity") +
    scale_fill_manual(values = c("darkblue", "royalblue", "lightblue", "tomato", "orange", "gold")) +
    scale_y_continuous("expected\nnb. of cells", limits = c(0, 125)) +
    .prettify("bw") + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        aspect.ratio = NULL,
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(1,0,0,2), "mm"))

#fns <- list.files("/users/helena/dropbox/portmac/kang", "ds10_ss[0-9];", full.names = TRUE)
res <- map(lapply(snakemake@input$res, readRDS), "tbl")
rmv <- vapply(res, is.null, logical(1))
res <- map(res[!rmv], mutate_if, is.factor, as.character) %>% 
    bind_rows %>% 
    mutate(E = (sim_mean.A + sim_mean.B) / 2) %>% 
    dplyr::filter(E > 0.1) %>% 
    setDT %>% split(by = c("sid", "i", "mid"), flatten = FALSE)

cd <- map_depth(res, 2, function(u) {sapply(u, nrow)
    truth <- setDF(select(u[[1]], c("sid", "i", "is_de")))
    pvals <- map(lapply(c("p_val", "p_adj.loc"), map, .x = u), setDF)
    dfs <- set_names(c(list(truth), pvals), c("truth", "pval", "padj"))
    do.call(COBRAData, dfs)
})

perf <- map_depth(cd, 2, 
    calculate_performance,
    binary_truth = "is_de",
    aspects = c("fdrtpr", "fdrtprcurve"))

gg_df <- map_depth(perf, 2, fdrtpr) %>% 
    map(bind_rows, .id = "i") %>% 
    bind_rows(.id = "sid") %>%
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>% 
    group_by(sid, thr, method) %>% 
    summarise_at(c("FDR", "TPR"), mean) %>% 
    mutate_at("method", factor, levels = names(.meth_cols))

p <- .plot_perf_points(gg_df, facet = "sid") +
    theme(plot.margin = unit(c(-1,0,0,2), "mm"))
p$facet$params$ncol <- nlevels(factor(gg_df$sid))

p <- cowplot::plot_grid(hists, p, 
    ncol = 1, axis = "lr", align = "v",
    rel_heights = c(1, 3))

ggsave(snakemake@output$fig, 
    width = 15, height = 7.1, units = "cm",
    dpi = 300, useDingbats = FALSE)

