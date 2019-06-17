library(RColorBrewer)
rs <- colorRampPalette(brewer.pal(9, "Reds")   [-seq_len(2)])(10)
gs <- colorRampPalette(brewer.pal(9, "Greens") [-seq_len(2)])(13)
ps <- colorRampPalette(brewer.pal(9, "Purples")[-seq_len(2)])(13)
method_colors <- c(
    "AD.logcounts" = "peru",
    "AD.vstcounts" = "sandybrown",
    "scDD.logcounts" = "maroon",
    "scDD.vstcounts" = "pink",
    "MAST.logcpm"       = rs[10],
    "MAST.logcpm_dr"    = rs[7],
    "MAST.logcounts"    = rs[4],
    "MAST.logcounts_dr" = rs[1],
    "MM-dream"    = "royalblue", 
    "MM-dream_dr" = "skyblue", 
    "MM-vst"      = "cyan",
    "MM-vst_dr"   = "aquamarine",
    "edgeR.sum(counts)"             = "black", 
    "limma-voom.sum(counts)"        = "grey40", 
    "limma-trend.sum(counts)"       = "grey80", 
    
    "limma-trend.median(cpm)"       = gs[10],
    "limma-trend.median(logcounts)" = gs[7],
    "limma-trend.median(vstcounts)" = gs[4],
    "limma-trend.median(scalecpm)"  = "gold",
    
    "limma-trend.sum(normcounts)"   = ps[13],
    "limma-trend.sum(logcounts)"    = ps[10],
    "limma-trend.sum(vstcounts)"    = ps[7],
    "limma-trend.sum(cpm)"          = ps[4],
    "limma-trend.sum(scalecpm)"     = ps[1],
    
    "limma-trend.mean(normcounts)"  = gs[13],
    "limma-trend.mean(logcounts)"   = gs[10],
    "limma-trend.mean(vstcounts)"   = gs[7],
    "limma-trend.mean(cpm)"         = gs[4],
    "limma-trend.mean(scalecpm)"    = gs[1])

# ------------------------------------------------------------------------------
# aesthetics for plotting
# ------------------------------------------------------------------------------
prettify <- function(theme = NULL, ...) {
    if (is.null(theme)) theme <- "classic"
    base <- paste0("theme_", theme)
    base <- getFromNamespace(base, "ggplot2")
    base(base_size = 8) + theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
        plot.title = element_text(face = "bold", hjust = 0),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(0.2, "cm"),
        strip.background = element_rect(fill = NA),
        panel.spacing = unit(0, "cm"),
        ...)
}

# ------------------------------------------------------------------------------
# plot TPR vs. FDR
# ------------------------------------------------------------------------------
.plot_perf_curves <- function(x) {
    suppressMessages(
        p <- plot_fdrtprcurve(x, pointsize = 1, linewidth = 0.6) + 
            scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.05, 0)) +
            scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.05, 0)) +
            prettify(theme = "bw", aspect.ratio = 1,
                axis.line = element_blank()) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))
    p$layers[[1]] <- NULL # remove vertical dashed lightgrey lines 
    p$layers[[1]]$aes_params$size <- 0.4 # dashed lines at thrs
    return(p)
}

.plot_perf_points <- function(df, 
    color_by = "method", colors = method_colors, facet = NULL) {
    p <- ggplot(df, aes_string(x = "FDR", y = "TPR", col = color_by)) +
        geom_vline(size = 0.2, lty = 2, aes(xintercept = thr)) + 
        geom_point(size = 1, alpha = 0.8) + 
        geom_line(size = 0.4, alpha = 0.4, show.legend = FALSE) +
        scale_color_manual(NULL, values = colors) +
        scale_x_sqrt(limits = c(0, 1), breaks = c(c(0.01, 0.1), seq(0, 1, 0.2)), 
            labels = function(x) format(x, drop0trailing = TRUE), expand = c(0, 0.05)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
        guides(col = guide_legend(ncol = 3,
            override.aes = list(size = 2, alpha = 1))) +
        prettify(theme = "bw", aspect.ratio = 1,
            legend.position = "bottom",
            legend.direction = "horizontal")
    if (!is.null(facet)) p <- p + facet_wrap(facet, 
        labeller = labeller(.multi_line = FALSE))
    return(p)
}

# ------------------------------------------------------------------------------
# helper to run dim. red. on max. of n cells per cluster-sample
#   x = a SingleCellExperiment
#   y = dimensionaly reduction method
#   n = max. nb. of cells to use per cluster-sample
# ------------------------------------------------------------------------------
.run_dr <- function(x, y = "TSNE", n = 100) {
    cells_use <- split(colnames(x), 
        list(x$cluster_id, x$sample_id)) %>% 
        map(function(u) sample(u, min(100, length(u))))
    fun <- paste0("run", y)
    fun <- getFromNamespace(fun, "scater")
    fun(x[, unlist(cells_use)])
}

# ------------------------------------------------------------------------------
# plot reduced dimensions from a list of SingleCellExperiment
# ------------------------------------------------------------------------------
.plot_dr <- function(x, y = NULL, color_by, lab) {
    if (is.null(y)) y <- reducedDimNames(x[[1]])[1]
    labs <- switch(y, TSNE = "t-SNE", y)
    labs <- paste(labs, "dim.", seq_len(2))
    get_labs <- function(n, lab) {
        labs <- paste0(n * 100, "% ", lab)
        labs[is.na(n)] <- "Reference data"
        return(labs)
    }
    df <- map(x, function(u) 
        data.frame(reducedDim(u, y), col = u[[color_by]])) %>% 
        bind_rows(.id = "facet") %>% 
        mutate_at("facet", function(u) {
            n <- as.numeric(u)
            labs <- get_labs(n, lab)
            lvls <- c("Reference data", get_labs(sort(unique(n)), lab))
            factor(labs, levels = lvls)
        })
    ggplot(df, aes(x = X1, y = X2, col = col)) + 
        geom_point(size = 0.1, alpha = 0.2) + 
        facet_wrap("facet", nrow = 1) +
        guides(color = guide_legend(color_by, 
            override.aes = list(alpha = 1, size = 2))) +
        labs(x = labs[1], y = labs[2]) + prettify() + 
        theme(aspect.ratio = 1, axis.line = element_blank())
}

# ------------------------------------------------------------------------------
# wrapper to extract method_ids, nb. of genes/cells, 
# and replicate number from snakemake output filenames
# ------------------------------------------------------------------------------
.parse_fns <- function(x) {
    list(
        method = factor(basename(gsub(";.*", "", x)), levels = names(method_colors)),
        n_genes = gsub(".*n_genes=([0-9]*).*", "\\1", x) %>% as.numeric,
        n_cells = gsub(".*n_cells=([0-9]*).*", "\\1", x) %>% as.numeric,
        n_samples = gsub(".*n_samples=([0-9]*).*", "\\1", x) %>% as.numeric,
        sim_rep = gsub(".*simrep=([0-9]+).*", "\\1", x) %>% factor,
        run_rep = gsub(".*runrep=([0-9]+).*", "\\1", x) %>% factor
    )
}

# ------------------------------------------------------------------------------
# wrapper to get metadata$gene_info slot from a SCE, assure correct formatting, 
# and add column is_de indicating whether a gene-cluster is DE (1) or not (0)
# ------------------------------------------------------------------------------
.get_gis <- function(fns) {
    rep <- gsub(".*simrep=([0-9]+).*", "\\1", fns)
    sce <- lapply(fns, readRDS)
    md <- map(sce, S4Vectors::metadata)
    gi <- map(md, function(u)
        u$gene_info %>% 
            mutate_if(is.factor, as.character) %>% select(-"logFC") %>% 
            mutate(sim_logFC = log2(sim_mean.B/sim_mean.A)) %>% 
            mutate(is_de = !category %in% c("ee", "ep")) %>% 
            mutate_at("is_de", as.integer))
    o <- order(as.numeric(rep))
    gi[o] %>% set_names(rep[o])
}

# ------------------------------------------------------------------------------
# wrapper to split results by sim_rep, run_rep & method,
# and return a list of tidy-format data.frames
# ------------------------------------------------------------------------------
.tidy_res <- function(fns, gis = NULL) {
    # load data
    bns <- basename(fns)
    res <- lapply(fns, readRDS)
    # extract results
    res <- map(res, "res")
    # remove failed runs
    rmv <- vapply(res, is.null, logical(1))
    bns <- bns[!rmv]
    res <- res[!rmv]
    # extract run-configuration from filenames
    pars <- lapply(.parse_fns(bns), function(u) 
        if (is.factor(u)) droplevels(u) else u)
    # construct data.frames
    cols_keep <- c("gene", "cluster_id", "p_val", "p_adj.loc", "p_adj.glb")
    if (!is.null(gis)) {
        dfs <- lapply(seq_along(res), function(i) {
            if ("logFC" %in% names(res[[i]])) 
                cols_keep <- c(cols_keep, "logFC")
            re <- res[[i]] %>% select(cols_keep)
            gi <- gis[[as.character(pars$sim_rep[i])]]
            df <- left_join(gi, re, c("gene", "cluster_id"))
            df[, names(pars)] <- map(pars, i)
            return(df)
        })
    } else {
        dfs <- lapply(seq_along(res), function(i) {
            if ("logFC" %in% names(res[[i]])) 
                cols_keep <- c(cols_keep, "logFC")
            df <- res[[i]] %>% select(cols_keep)
            df[, names(pars)] <- map(pars, i)
            return(df)
        })
    }
    dfs %>% bind_rows %>% 
        dplyr::rename(est_logFC = logFC) %>% 
        mutate_at(c("sim_rep", "run_rep"), function(u)
            factor(as.character(u), levels = sort(as.numeric(levels(u)))))
}

# ------------------------------------------------------------------------------
# removes lowest contribution row/column until all entries >= n, 
# or until dimensions match input argument 
# (used to filter clusters/samples to use for simulation)
# ------------------------------------------------------------------------------
.filter_matrix <- function(m, n = 100, dim = NULL) {
    if (!is.null(dim)) {
        test <- expression(nrow(m) > dim[1] | ncol(m) > dim[2])
        test_row <- expression(nrow(m) <= dim[1])
        test_col <- expression(ncol(m) <= dim[2])
    } else {
        test <- expression(any(m < n))
        test_row <- expression(all(m[rmv, ] >= n))
        test_col <- expression(all(m[, rmv] >= n))
    }
    while (eval(test)) {
        s <- sum(m)
        rows <- rowSums(m) / s
        cols <- colSums(m) / s
        x <- c(rows, cols)
        rmv <- names(which.min(x))
        y <- TRUE
        while (y) {
            if (rmv %in% rownames(m)) {
                if (eval(test_row)) {
                    x <- x[names(x) != rmv]
                    rmv <- names(which.min(x)) 
                } else {
                    m <-  m[rownames(m) != rmv, , drop = FALSE]
                    y <- FALSE
                }
            } else {
                if (eval(test_col)) {
                    x <- x[names(x) != rmv]
                    rmv <- names(which.min(x))
                } else {
                    m <- m[, colnames(m) != rmv, drop = FALSE]
                    y <- FALSE
                }
            }
        }
    }
    return(m)
}

# ------------------------------------------------------------------------------
# filter a SCE provided with character vectors of clusters & samples to keep
# ------------------------------------------------------------------------------
.filter_sce <- function(x, y) {
    cells_keep <- data.frame(cell = colnames(x), colData(x)) %>% 
        dplyr::filter(cluster_id %in% y[[1]] & sample_id %in% y[[2]]) %>% 
        select("cell") %>% unlist %>% as.character
    x <- x[, cells_keep]
    x$cluster_id <- droplevels(x$cluster_id)
    x$sample_id <- droplevels(x$sample_id)
    ei <- metadata(x)$experiment_info
    ei <- ei[ei$sample_id %in% levels(x$sample_id), ]
    ei$sample_id <- droplevels(ei$sample_id)
    ei$group_id <- droplevels(ei$group_id)
    metadata(x)$experiment_info <- ei
    return(x)
}

# ------------------------------------------------------------------------------
# update SCE colData such that only levels of existent factors are retained
# ------------------------------------------------------------------------------
.update_cd <- function(sce) {
    colData(sce) <- as.data.frame(colData(sce)) %>% 
        mutate_if(is.factor, droplevels) %>% 
        DataFrame(row.names = colnames(sce))
    return(sce)
}

# ------------------------------------------------------------------------------
# calculate performance using iCOBRA
#   - df: data.frame containing p_val & p_adj.X
#   - gi: data.frame as returned by simData() 
#         in metadata()$gene_info slot;
#         containing column is_de
# ------------------------------------------------------------------------------
.calc_perf <- function(x, facet = "none", ...) {
    library(iCOBRA)
    cd <- COBRAData(pval = x$p_val, padj = x$p_adj, truth = x$truth)
    suppressMessages(calculate_performance(...,
        cd, binary_truth = "is_de", 
        splv = facet, maxsplit = Inf, 
        aspects = c("fdrtpr", "fdrtprcurve"))) %>% 
        prepare_data_for_plot(facetted = (facet != "none"))
}

# ------------------------------------------------------------------------------
# for ea. k in 1:maxrank, count the nb. of genes occurring ea. nb. of times
# > 3 column data.frame: nbr_occ | nbr_genes | k 
#   (for a given row, among the top-k genes from each column of mtx, 
#    nbr_genes occur exactly nbr_occ times)
# ------------------------------------------------------------------------------
.calc_occs <- function(mtx, maxrank) {
    maxrank <- min(maxrank, nrow(mtx))
    ks <- seq_len(maxrank)
    if (ncol(mtx) > 1) {
        M <- matrix(0, max(mtx[ks, ]), maxrank)
        for (i in seq_len(ncol(mtx)))
            M[cbind(mtx[ks, i], ks)] <- M[cbind(mtx[ks, i], ks)] + 1
        M <- M[rowSums(M) != 0, ]
        M <- t(apply(M, 1, cumsum))
        M2 <- matrix(0, ncol(mtx), maxrank)
        for (i in 1:nrow(M))
            M2[cbind(M[i, ], seq_len(ncol(M)))] <- 
                M2[cbind(M[i, ], seq_len(ncol(M)))] + 1
        M2 <- M2 %>% reshape2::melt() %>%
            dplyr::rename(k = Var2, nbr_genes = value, nbr_occ = Var1)
        M2 %>% dplyr::arrange(k, nbr_occ) %>%
            dplyr::mutate(nbr_cols = ncol(mtx))
    } else {
        NULL
    }
}

# ------------------------------------------------------------------------------
# calculate partial (cumulative) AUCs. 
# (assumes that x variable = k, y variable = nbr_genes
# ------------------------------------------------------------------------------
.calc_aucs <- function(x) {
    x %>% dplyr::mutate(dx = c(k[1], diff(k)),
        dy = c(nbr_genes[1], diff(nbr_genes)),
        ys = c(0, nbr_genes[-length(nbr_genes)])) %>%
        dplyr::mutate(AUC = cumsum(dx * dy/2 + dx * ys)) %>%
        dplyr::mutate(AUCs = AUC/(k^2/2))
}

# ------------------------------------------------------------------------------
# get all subclusters from an hclust object
# ------------------------------------------------------------------------------
.get_subcl <- function(hcl) {
    m <- hcl$merge
    labs <- hcl$labels
    L <- list()
    for (i in seq_len(nrow(m))) {
        tmp <- c()
        if (m[i, 1] < 0) tmp <- c(tmp, labs[-m[i, 1]])
        else tmp <- c(tmp, L[[m[i, 1]]])
        if (m[i, 2] < 0) tmp <- c(tmp, labs[-m[i, 2]])
        else tmp <- c(tmp, L[[m[i, 2]]])
        L[[i]] <- sort(tmp)
    }
    L
}
