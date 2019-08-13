# load packages
suppressMessages({
    library(circlize)
})

# load data & results
sce <- readRDS(file.path("MAGL", "output", "MAGL-SCE.rds"))

.deg <- function(res, fdr = 0.05, lfc = 1)
    map(map(res$table[[1]], dplyr::filter, p_adj.loc < fdr, abs(logFC) > 1), pull, "gene") 

# construct data.frame for plotting
ids <- c("cluster_id", "group_id", "sample_id")
cd_df <- data.frame(colData(sce), do.call(cbind, reducedDims(sce)))
cd_df <- cd_df[sample(nrow(cd_df), nrow(cd_df)), ] # randomize rows

# store IDs & nb. of levels
nk <- length(kids <- purrr::set_names(levels(sce$cluster_id)))
ns <- length(sids <- purrr::set_names(levels(sce$sample_id)))
ng <- length(gids <- purrr::set_names(levels(sce$group_id)))

# specify color palettes
cluster_id_pal <- muscat:::.cluster_colors[c(1, 3, 5, 7, 13, 17, 19, 11, 1)]
sample_id_pal <- CATALYST:::.cluster_cols[seq_len(ns) + nk]
group_id_pal <- c("royalblue", "orange")

# round limits of x to nearest n
.get_brks <- function(x, n) {
    n <- 1 / n
    min <- floor(min(x) * n) / n
    max <- ceiling(max(x) * n) / n
    c(0, min, max)
}

# color palettes for heatmap, cluster, sample, group IDs, and # cells
hm_cols <- c("grey95", "blue3", "red3")
kcols <- set_names(CATALYST:::.cluster_cols[seq_len(nk)], kids)
scols <- set_names(CATALYST:::.cluster_cols[seq_len(ns) + nk], sids)
gcols <- set_names(c("royalblue", "orange"), gids)

n_cells <- c(t(metadata(pb)$n_cells))
max <- .get_brks(n_cells, 100)[3]
ncols <- colorRamp2(c(0, max), c("white", "black"))

lgd_aes <- list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6))

# column annoation 
col_df <- data.frame(n_cells, 
    cluster_id = rep(kids, each = ns), 
    group_id = rep(pb$group_id, nk))
col_cols <- list(
    n_cells = ncols, cluster_id = kcols, 
    sample_id = scols, group_id = gcols)
col_anno <- HeatmapAnnotation(
    df = col_df, col = col_cols,
    show_annotation_name = FALSE,
    annotation_legend_param = lgd_aes)

# heatmap wrapper
.plot_hm <- function(sce, col, brks, col_anno = NULL, row_anno = NULL, 
    row_split = NULL, row_nms = FALSE, col_title = NULL, ...)
    Heatmap(assay(sce), col, name = "logFC", 
        column_title = col_title, cluster_columns = FALSE, 
        show_row_names = row_nms, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = TRUE,
        row_split = row_split, column_split = sce$cluster_id,
        top_annotation = col_anno, left_annotation = row_anno, 
        heatmap_legend_param = c(lgd_aes, list(at = brks)),
        row_names_side = "left", ...)
