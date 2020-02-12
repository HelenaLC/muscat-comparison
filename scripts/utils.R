.typ_cols <- c(pb = "grey50", mm = "violet", ad = "lightgreen", mast = "orange", scdd = "steelblue")
.typ_labs <- c(pb = "PB", mm = "MM", ad = "AD", mast = "MAST", scdd = "scDD")

.cat_cols <- c("royalblue", "cornflowerblue", "red3", "tomato", "orange", "gold")
names(.cat_cols) <- c("ee", "ep", "de", "dp", "dm", "db")

.meth_cols <- c(
    "edgeR.sum.counts" = "#000000",
    "edgeR.sum.scalecpm" = "#C6C6C6",
    
    "scDD.logcounts"    = "#0056B2",
    "scDD.vstresiduals" = "#009EF6", 
    
    "limma-voom.sum.counts" = "#43A047",
    "limma-trend.mean.logcounts"    = "#45BF55",
    "limma-trend.mean.vstresiduals" = "#B5E655",
    
    "MAST.logcounts"    = "#FFC300",
    
    "MM-dream"  = "#005E5C",
    "MM-dream2" = "#00ABA6",
    "MM-nbinom" = "#00E3DD",
    "MM-vst"    = "#95EEE8",
    
    "AD-gid.logcounts"    = "#9A41B3",
    "AD-gid.vstresiduals" = "#FFA9FF",
    "AD-sid.logcounts"    = "#E56D4B",
    "AD-sid.vstresiduals" = "#FBB6A2")

#cols <- .meth_cols
#hist(seq_along(cols), breaks = c(seq_along(cols) - 0.5, length(cols) + 0.5), col = cols)

.read_res <- function(fns, slot = "tbl") {
    res <- map(lapply(fns, readRDS), slot)
    rmv <- vapply(res, function(u) 
        is.null(u) | inherits(u, "error"), 
        logical(1))
    res <- res[!rmv]
    if (slot == "tbl")
        res <- map(res, mutate_if, is.factor, as.character) %>% 
            bind_rows %>% mutate_if(is.character, as.factor) %>% 
            mutate_at("category", factor, levels = muscat:::cats) %>% 
            mutate_at("mid", factor, levels = names(.meth_cols)) %>% 
            mutate_if(is.factor, droplevels)
    return(res)
}

.filter_matrix <- function(m, n = 100) {
    while (any(m < n)) {
        # get candidate rows/cols for removal
        i <- m < n
        r <- apply(i, 1, any)
        c <- apply(i, 2, any)
        # get smallest row/col
        rs <- rowSums(m)
        cs <- colSums(m)
        r <- which(r)[which.min(rs[r])]
        c <- which(c)[which.min(cs[c])]
        # priorities removal of rows over cols
        if (rs[r] <= cs[c]) {
            m <- m[-r, , drop = FALSE]
        } else {
            m <- m[, -c, drop = FALSE]
        }
        if (any(dim(m) == 1)) 
            break
    }
    return(m)
}

.update_sce <- function(sce) {
    # update colData
    cd <- as.data.frame(colData(sce))
    cd <- mutate_if(cd, is.factor, droplevels) 
    colData(sce) <- DataFrame(cd, row.names = colnames(sce))
    # update metadata
    ei <- metadata(sce)$experiment_info
    ei <- ei[ei$sample_id %in% levels(sce$sample_id), ]
    ei <- mutate_if(ei, is.factor, droplevels)
    metadata(sce)$experiment_info <- ei
    return(sce)
}

.filter_sce <- function(sce, kids, sids) {
    cs1 <- sce$cluster_id %in% kids
    cs2 <- sce$sample_id %in% sids
    sce <- sce[, cs1 & cs2]
    sce <- .update_sce(sce)
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
        legend.key.size = unit(2, "mm"),
        strip.background = element_rect(fill = NA),
        plot.margin = unit(rep(1, 4), "mm"),
        panel.spacing = unit(0, "mm"),
        legend.margin = margin(0,0,1,0,"mm"),
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
            guides(col = guide_legend(ncol = 4,
                override.aes = list(size = 2, alpha = 1))) +
            .prettify(theme = "bw", legend.position = "bottom"))

