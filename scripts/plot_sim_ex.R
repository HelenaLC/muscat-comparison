source(snakemake@config$utils)

suppressMessages({
    library(ComplexHeatmap)
    library(magrittr)
    library(muscat)
    library(purrr)
    library(SingleCellExperiment)
})

set.seed(3004)

sce <- readRDS("data/raw_data/kang.rds")
nk <- nlevels(sce$cluster_id)
ns <- nlevels(sce$sample_id)

sim <- simData(sce, nrow(sce), 2*nk*ns*100,
    p_dd = c(0.7, 0.1, rep(0.05, 4))) %>% 
    scran::computeSumFactors %>% 
    normalize

m <- match(colnames(mat), sim$sample_id)
groups <- sim$group_id[m]

ms <- aggregateData(sim, "logcounts", fun = "mean") %>% 
    assays %>% as.list %>% 
    map(data.frame, gene = rownames(sim)) %>%  
    bind_rows(.id = "cluster_id") %>%
    set_rownames(sprintf("%s(%s)", .$gene, .$cluster_id)) %>% 
    select(-c("gene", "cluster_id"))

gi <- metadata(sim)$gene_info %>% 
    mutate(id = sprintf("%s(%s)", gene, cluster_id)) %>% 
    mutate(logFC = log2(sim_mean.A/sim_mean.B)) %>% 
    mutate_at("category", plyr::revalue, c("de" = "ds")) %>% 
    mutate_at("category", factor, labels = toupper(levels(.$category)))

gs <- sample_n(group_by(gi, category), 8)

mat <- as.matrix(ms[gs$id, ])
mat <- muscat:::.z_norm(mat)

hm_cols <- c("darkgrey", "white", "slateblue")#c("skyblue", "cornflowerblue", "royalblue", "black", "orange3", "orange", "gold")
cat_cols <- setNames(c("blue3", "cornflowerblue", "red3", "tomato", "orange", "gold"),levels(gi$category))
lfc_cols <- circlize::colorRamp2(c(-4,0,4), c("limegreen", "white", "tomato"))
grp_cols <- setNames(scales::hue_pal()(2), levels(sim$group_id))

row_anno <- rowAnnotation(
    df = data.frame(
        logFC = gs$logFC,
        category = gs$category),
    col = list(
        logFC = lfc_cols,
        category = cat_cols),
    gp = gpar(col = "white"),
    simple_anno_size = unit(4, "mm"),
    show_annotation_name = FALSE)

col_anno <- columnAnnotation(
    df = data.frame(
        group_id = groups),
    col = list(
        group_id = grp_cols),
    gp = gpar(col = "white"),
    simple_anno_size = unit(4, "mm"),
    show_annotation_name = FALSE)

hm <- Heatmap(mat, 
    col = hm_cols,
    name = "z-normazlied\nmean expr.",
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    rect_gp = gpar(col = "white"), 
    left_annotation = row_anno,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_split = gs$category,
    column_split = groups,
    row_title = NULL,
    column_title = NULL,
    top_annotation = col_anno)
   
hm <- grid::grid.grabExpr(draw(hm,
    heatmap_row_names_gp = gpar(fontsize = 4),
    heatmap_column_names_gp = gpar(fontsize = 8),
    legend_title_gp = gpar(fontsize = 8),
    legend_labels_gp = gpar(fontsize = 6)))

ggsave("/Users/helena/Desktop/hm.pdf", hm,
    width = 15, height = 14, units = "cm",
    dpi = 300, useDingbats = FALSE)

