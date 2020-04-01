config <- file.path("..", "config.yaml")
config <- yaml::read_yaml(config)
source(file.path("..", config$utils))

# load packages
suppressMessages({
    library(circlize)
    library(ComplexHeatmap)
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(RColorBrewer)
    library(SingleCellExperiment)
})

# get experiment metadata table
ei <- metadata(sce)$experiment_info

# construct data.frame of cell metadata & reduced dimensions
ids <- c("cluster_id", "group_id", "sample_id")
cd_df <- data.frame(cell_id = seq_len(ncol(sce)),
    colData(sce), do.call(cbind, reducedDims(sce)))
cd_df <- cd_df[sample(nrow(cd_df), nrow(cd_df)), ] # randomize rows

# store IDs & nb. of levels
nk <- length(names(kids) <- kids <- levels(sce$cluster_id))
ns <- length(names(sids) <- sids <- levels(sce$sample_id))
ng <- length(names(gids) <- gids <- levels(sce$group_id))

# color palettes for cluster, sample and group IDs -----------------------------
pal <- muscat:::.cluster_colors
cluster_id_pal <- c(
    "Inhib. Neuron" = pal[1],
    "Excit. Neuron" = pal[11],
    "Astrocytes" = pal[3],
    "Endothelial" = pal[5],
    "OPC" = pal[14],
    "Oligodendrocytes" = pal[13],
    "Microglia" = pal[7],
    "CPE cells" = pal[20])
sample_id_pal <- c(
    brewer.pal(9, "Blues")[c(3,5,7,9)],
    brewer.pal(9, "Oranges")[c(3,5,7,9)])
names(sample_id_pal) <- ei$sample_id[c(
    which(ei$group_id == "Vehicle"),
    which(ei$group_id == "LPS"))]
group_id_pal <- c(Vehicle = "royalblue", LPS = "orange")

# save color legends to .pdf ---------------------------------------------------
for (id in ids) {
    fn <- sprintf(file.path("figures", "lgd_%s.pdf"), id)
    pal <- get(paste0(id, "_pal"))
    foo <- ggplot(cd_df, aes_string("ident", "ident", col = id)) + 
        scale_color_manual(values = pal) + geom_point() + .prettify()
    ggsave(fn, get_legend(foo), 
        width = 3, height = 4, units = "cm", 
        dpi = 300, useDingbats = FALSE)
}

# get axis limits by rounding limits of x to nearest n -------------------------
.get_brks <- function(x, n) {
    n <- 1 / n
    min <- floor(min(x) * n) / n
    max <- ceiling(max(x) * n) / n
    c(0, min, max)
}

# wrapper to plot prettified dimension reduction plot --------------------------
.plot_dr <- function(df, col, pal)
    ggplot(df, aes_string(x = "UMAP_1", y = "UMAP_2", col = col)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    geom_point_rast(size = 0.1, alpha = 0.1) + 
    scale_color_manual(values = pal) + 
    theme_void() + theme(aspect.ratio = 1)

# cluster-wise logFC heatmap wrapper -------------------------------------------
.hm_cols <- c("navy", "cyan", "violet")
.hm_lgd_aes <- list(
    title_gp = gpar(fontsize = 8), 
    labels_gp = gpar(fontsize = 6))
.plot_hm <- function(sce, col, brks, col_anno = NULL, row_anno = NULL, 
    row_split = NULL, row_nms = FALSE, col_title = NULL, ...)
    Heatmap(assay(sce), col, name = "logFC", 
        column_title = col_title, cluster_columns = FALSE, 
        show_row_names = row_nms, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = TRUE,
        row_split = row_split, column_split = sce$cluster_id,
        top_annotation = col_anno, left_annotation = row_anno, 
        heatmap_legend_param = c(.hm_lgd_aes, list(at = brks)),
        row_names_side = "left", ...)

# wrapper to extract DS genes per cluster --------------------------------------
.get_ds_gs <- function(u, fdr = 0.05, lfc = 1, simplify = FALSE) {
    u <- map(u, filter, p_adj.loc < fdr)
    if ("logFC" %in% colnames(u))
        u <- filer(u, abs(logFC) > lfc)
    if (simplify) unique(unlist(map(u, "gene"))) else u
}
    