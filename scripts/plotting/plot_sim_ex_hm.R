# load packages
suppressMessages(
    suppressPackageStartupMessages({
        library(circlize)
        library(ComplexHeatmap)
        library(dplyr)
        library(ggplot2)
        library(purrr)
        library(muscat)
        library(scales)
        library(scater)
        library(SingleCellExperiment)
        library(viridis)
    })
)

# load reference dataset
sce <- readRDS(file.path("data", "raw_data", "kang-sim.rds"))
nk <- nlevels(sce$cluster_id)
ns <- nlevels(sce$sample_id)

# simulate data
set.seed(1)
sim <- simData(sce, 
    n_genes = 5e3,
    n_cells = 2*100*nk*ns,
    p_dd = c(0.75, rep(0.05, 5)))

# calculate means by cluster-sample
sim <- normalize(sim)
pb <- aggregateData(sim, 
    by = c("cluster_id", "sample_id"), 
    assay = "logcounts", fun = "mean")

# sample n genes per category
n <- 10
gi <- metadata(sim)$gene_info %>% 
    dplyr::mutate(logFC = log2(sim_mean.B/sim_mean.A))
gi_by_cat <- split(gi, gi$category)
cats <- names(gi_by_cat)
gs <- lapply(cats, function(c) {
    gs <- gi_by_cat[[c]]
    gs[sample(seq_len(nrow(gs)), n), c("gene", "cluster_id")]
}) %>% set_names(cats) %>% 
    bind_rows(.id = "category")

# prep. matrix for plotting
cells_by_ks <- muscat:::.split_cells(sim)
ms <- sapply(seq_len(nrow(gs)), function(i) {
    g <- gs[i, "gene"]
    k <- gs[i, "cluster_id"]
    assays(pb)[[k]][g, ]
}) %>% t %>% set_rownames(with(gs, sprintf("%s(%s)", gene, cluster_id)))

# scale expression means
ms0 <- muscat:::.scale(ms)

# aesthetics
lgd_aes <- list(
    title_gp = gpar(fontsize = 8), 
    labels_gp = gpar(fontsize = 6))

# 1st row annotation: simulated logFC
lfcs <- left_join(gs, gi, by = c("gene", "cluster_id")) %>% pull("logFC")
min <- floor(min(lfcs) * 10) / 10
max <- ceiling(max(lfcs) * 10) / 10
lfc_anno <- rowAnnotation(
    df = data.frame(logFC = lfcs),
    col = list(logFC = colorRamp2(c(min, 0, max), 
        c("tomato", "white", "limegreen"))),
    gp = gpar(col = "white"),
    annotation_legend_param = lgd_aes)

# 2nd row annotation: gene category
cat_anno <- rowAnnotation(
    df = data.frame(category = gs$category),
    col = list(category = c("ee" = "blue3", "ep" = "cornflowerblue",
        "de" = "red3", "dp" = "tomato", "dm" = "orange", "db" = "gold")),
    gp = gpar(col = "white"),
    annotation_legend_param = lgd_aes)

# column annotation: group_id
ei <- metadata(sim)$experiment_info
m <- match(colnames(ms), ei$sample_id)
col_anno <- columnAnnotation(
    df = data.frame(group_id = ei$group_id[m]),
    col = list(group_id = setNames(hue_pal()(2), c("A", "B"))),
    gp = gpar(col = "white"),
    annotation_legend_param = lgd_aes)

# plot heatmap of expression means by cluster-sample
hm <- Heatmap(ms0, 
    col = viridis(10, option = "magma"),
    #col = rev(brewer.pal(11, "PuOr")),
    name = "scaled mean\nexpression",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.1),
    top_annotation = col_anno,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = lgd_aes)

# combine panels
p <- grid.grabExpr(draw(lfc_anno + cat_anno + hm))
ggsave("figures/sim_ex_hm.pdf", p,
    width = 15, height = 15, units = "cm",
    dpi = 300, useDingbats = FALSE)

