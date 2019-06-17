# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(cowplot)
        library(dplyr)
        library(ggplot2)
        library(ggtree)
        library(tidytree)
    }))

# source utils
source(snakemake@input$utils)

# load results
res <- readRDS(snakemake@input$res)

ggt <- ggtree(as.phylo(res$tree_avg), branch.length = "none")
for (t in seq_along(res$subc_avg)) {
    tryCatch({
        i <- ggtree::MRCA(ggt, res$subc_avg[[t]])
        if (res$stab_sco[t] > 0)
            ggt$data[i, "label"] <- round(res$stab_sco[t], 2)
    }, error = function(e) NULL, warning = function(w) NULL)
}

# add node labels & prettify
ggt <- ggt + 
    geom_label2(aes(subset = !isTip, label = label), size = 2) +
    geom_tiplab(aes(angle = 90), size = 2, hjust = 1) +
    scale_x_reverse() + coord_flip() + xlim_tree(16) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

tipo <- ggt$data %>% filter(isTip) %>% arrange(y)
anno <- readRDS("metadata/anno_tbl.rds")
anno <- anno[match(tipo$label, anno$id), ] %>% 
    mutate_at("inp", function(u) 
        factor(gsub("^log", "", u),
            levels = c("counts", "cpm", "normcounts", "vstcounts")))

cols <- list(
    typ = setNames(brewer.pal(3, "Reds"),
        paste0(c("group", "sample", "cell"), "-level")),
    inp = setNames(brewer.pal(nlevels(anno$inp), "Blues"), levels(anno$inp)),
    log = setNames(c("lightgreen", "lightcoral"), c("TRUE", "FALSE")),
    agg = setNames(
        c(rev(brewer.pal(3, "Purples")), "lightgrey"),
        c("sum", "mean", "median", "none")))

# generate annotation heatmaps
hms <- lapply(seq_along(cols), function(i) {
    ggplot(anno) + geom_tile(aes_string(x = seq_len(nrow(anno)), y = 1, fill = names(cols)[i])) + 
        geom_vline(xintercept = (seq_len(nrow(anno)-1)) + 0.5, 
            linetype = "solid", color = "white", size = 0.3) + 
        scale_fill_manual(values = cols[[i]]) + 
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + prettify() +
        theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(0, 3, 0, 3), "mm"))
})

# arrange panels
lgd <- lapply(hms, get_legend)
hms <- lapply(hms, function(p) p + theme(legend.position = "none"))
p1 <- plot_grid(
    plotlist = c(list(ggt), hms),
    ncol = 1, rel_heights = c(8, rep(0.5, length(hms))))
p2 <- plot_grid(plotlist = lgd, ncol = 1, axis = "l", align = "v")
p <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 0.2))

ggsave(snakemake@output$fig, p,
    width = 15, height = 10, units = "cm",
    dpi = 300, useDingbats = FALSE)

