.meth_cols <- c(
    "edgeR.sum(counts)" = "#000000",
    "edgeR.sum(scalecpm)" = "#C6C6C6",
    "limma-voom.sum(counts)" = "#5E5E5E",
    "limma-trend.mean(logcounts)"    = "#007700",
    "limma-trend.mean(vstresiduals)" = "#00B161",
    
    "MM-dream"  = "#336A7C",
    "MM-nbinom" = "#00E3DD",
    "MM-vst"    = "#3BB3B7",
    "scDD.logcounts"    = "royalblue",
    "scDD.vstresiduals" = "cornflowerblue", 
    
    "AD-gid.logcounts"    = "#B03060",
    "AD-gid.vstresiduals" = "#FFABD3",
    "AD-sid.logcounts"    = "#9E6A92",
    "AD-sid.vstresiduals" = "#E29CF2",
    "MAST.logcounts"    = "#A06A50")

#cols <- .meth_cols
#hist(seq_along(cols), breaks = c(seq_along(cols) - 0.5, length(cols) + 0.5), col = cols)

.filter_sce <- function(sce, kids, sids) {
    cs1 <- sce$cluster_id %in% kids
    cs2 <- sce$sample_id %in% sids
    sce <- sce[, cs1 & cs2]
    sce$cluster_id <- droplevels(sce$cluster_id)
    sce$sample_id <- droplevels(sce$sample_id)
    ei <- sce@metadata$experiment_info
    ei <- ei[ei$sample_id %in% levels(sce$sample_id), ]
    ei$sample_id <- droplevels(ei$sample_id)
    ei$group_id <- droplevels(ei$group_id)
    sce@metadata$experiment_info <- ei
    return(sce)
}

.update_cd <- function(sce) {
    colData(sce) <- as.data.frame(colData(sce)) %>% 
        mutate_if(is.factor, droplevels) %>% 
        DataFrame(row.names = colnames(sce))
    return(sce)
}

.prettify <- function(theme = NULL, ...) {
    if (is.null(theme)) theme <- "classic"
    base <- paste0("theme_", theme)
    base <- getFromNamespace(base, "ggplot2")
    base(base_size = 8) + theme(
        aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
        plot.title = element_text(face = "bold", hjust = 0),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(0.2, "cm"),
        strip.background = element_rect(fill = NA),
        panel.spacing = unit(0, "cm"),
        ...)}

.plot_perf_points <- function(df, color_by = "method", facet = "splitval")
    suppressMessages(
        ggplot(df, aes_string(x = "FDR", y = "TPR", col = color_by)) +
            facet_wrap(facet, labeller = labeller(.multi_line = FALSE)) +
            geom_vline(size = 0.2, lty = 2, aes(xintercept = thr)) + 
            geom_point(size = 1, alpha = 0.8) + 
            geom_line(size = 0.4, alpha = 0.4, show.legend = FALSE) +
            scale_color_manual(NULL, values = .meth_cols) +
            scale_x_sqrt(limits = c(0, 1), breaks = c(c(0.01, 0.1), seq(0.2, 1, 0.2)), 
                labels = function(x) format(x, drop0trailing = TRUE), expand = c(0, 0.05)) +
            scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
            guides(col = guide_legend(ncol = 3,
                override.aes = list(size = 2, alpha = 1))) +
            .prettify(theme = "bw", legend.position = "bottom"))

