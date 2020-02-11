config <- yaml::read_yaml("config.yaml")
source(config$utils)
source(file.path("MAGL", "code", "utils.R"))
ids <- c("cluster_id", "group_id", "sample_id")
dir <- "/users/helena/dropbox/phd/manuscript/fig4"

# load packages
suppressMessages({
    library(AnnotationDbi)
    library(circlize)
    library(ComplexHeatmap)
    library(dplyr)
    library(edgeR)
    library(ggplot2)
    library(M3C)
    library(Matrix)
    library(msigdbr)
    library(muscat)
    library(RColorBrewer)
    library(scran)
    library(SingleCellExperiment)
    library(topGO)
    library(org.Mm.eg.db)
    library(purrr)
    library(UpSetR)
})

# load data & results
sce <- readRDS(file.path("MAGL", "output", "MAGL-SCE.rds"))
sce$group_id <- factor(sce$group_id, 
    levels = levels(sce$group_id), 
    labels = c("vehicle", "LPS"))

pb <- aggregateData(sce)
res <- pbDS(pb)
deg_by_k <- .deg(res, fdr = 0.05, lfc = 1)
deg <- unique(unlist(deg_by_k))

# color legends ----------------------------------------------------------------
for (id in ids) {
    fn <- sprintf(file.path(dir, "magl_lgd-%s.pdf"), id)
    pal <- get(paste0(id, "_pal"))
    foo <- ggplot(cd_df, aes_string("ident", "ident", col = id)) + 
        scale_color_manual(values = pal) + geom_point() + .prettify()
    ggsave(fn, get_legend(foo), 
        width = 3, height = 4, units = "cm", 
        dpi = 300, useDingbats = FALSE)
}

# upset plot of DE genes -------------------------------------------------------
pdf(file.path(dir, "magl-upset.pdf"), width = 15/2.54, height = 9/2.54, onefile = FALSE)
upset(fromList(deg_by_k), sets = rev(kids)[-c(1, 2)],
    mb.ratio = c(0.6, 0.4), nintersects = 30,  
    text.scale = 0.8, point.size = 2, line.size = 0.4)
dev.off()

# relative cluster abundances --------------------------------------------------
# p1 <- table(sce$cluster_id, sce$sample_id) %>% 
#     prop.table(2) %>% unclass %>% melt %>% 
#     set_colnames(c("cluster_id", "sample_id", "frequency")) %>% 
#     dplyr::mutate(group_id = sce$group_id[match(.$sample_id, sce$sample_id)]) %>% 
#     ggplot(aes(x = sample_id, y = frequency, fill = cluster_id)) +
#     geom_bar(stat = "identity", width = 1, size = 0.2, col = "white") + 
#     facet_wrap("group_id", ncol = 1, scales = "free_y") +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
#     scale_fill_manual(values = cluster_id_pal) +
#     coord_flip() +
#     .prettify("bw") + theme(aspect.ratio = NULL, panel.spacing = unit(1, "mm"), 
#         strip.background = element_blank(), strip.text = element_blank(),
#         panel.grid = element_blank(), plot.margin = margin(1,3,0,0,"mm"))
# 
# sub <- sce[, !sce$cluster_id %in% c("Excit. Neuron", "Inhib. Neuron")]
# p2 <- table(sub$cluster_id, sub$sample_id) %>% 
#     prop.table(2) %>% unclass %>% melt %>% 
#     set_colnames(c("cluster_id", "sample_id", "frequency")) %>% 
#     dplyr::mutate(group_id = sce$group_id[match(.$sample_id, sce$sample_id)]) %>% 
#     ggplot(aes(x = sample_id, y = frequency, fill = cluster_id)) +
#     geom_bar(stat = "identity", width = 1, size = 0.2, col = "white") + 
#     facet_wrap("group_id", ncol = 1, scales = "free_y") +
#     scale_x_discrete(NULL, expand = c(0, 0)) +
#     scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
#     scale_fill_manual(values = cluster_id_pal) +
#     coord_flip() +
#     .prettify("bw") + theme(aspect.ratio = NULL, panel.spacing = unit(1, "mm"), 
#         strip.background = element_blank(), strip.text = element_blank(),
#         axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
#         panel.grid = element_blank(), plot.margin = margin(1,3,0,0,"mm"))
# 
# p3 <- table(sce$sample_id) %>% 
#     unclass %>% melt %>% 
#     set_colnames(c("sample_id", "n_cells")) %>% 
#     dplyr::mutate(group_id = sce$group_id[match(.$sample_id, sce$sample_id)]) %>% 
#     ggplot(aes(x = sample_id, y = n_cells, fill = group_id)) +
#     geom_bar(stat = "identity", width = 1, col = "white", size = 0.2) + 
#     facet_wrap("group_id", ncol = 1, scales = "free_y") +
#     scale_x_discrete(NULL, expand = c(0, 0)) +
#     scale_y_continuous("nb. of cells", limits = c(0, 5e3), breaks = seq(0, 5e3, 2500), expand = c(0,0)) +
#     scale_fill_manual(values = group_id_pal) +
#     coord_flip() +
#     .prettify("bw") + theme(aspect.ratio = NULL, panel.spacing = unit(1, "mm"), 
#         strip.background = element_blank(), strip.text = element_blank(),
#         axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
#         plot.margin = margin(1,4,0,0,"mm"))
# 
# p <- cowplot::plot_grid(nrow = 1, 
#     p1 + theme(legend.position = "none"), 
#     p2 + theme(legend.position = "none"), 
#     p3 + theme(legend.position = "none"),
#     plot_grid(get_legend(p1), get_legend(p3), ggplot() + theme_nothing(), 
#         ncol = 1, align = "v", axis = "l", rel_heights = c(3,1.5,2.5)),
#     align = "h", axis = "tl",
#     rel_widths = c(3,2,2,1.2))
# 
# ggsave(file.path(dir, "magl-cluster_freqs.pdf"), p,
#     width = 15, height = 5, units = "cm",
#     dpi = 300, useDingbats = FALSE) 

# UMAP colored by cluster, sample & group ID -----------------------------------
.plot_dr <- function(df, col, pal)
    ggplot(df, aes_string(x = "UMAP_1", y = "UMAP_2", col = col)) +
    geom_point_rast(size = 0.1, alpha = 0.1, show.legend = FALSE) + 
    scale_color_manual(values = pal) + 
    theme_void() + theme(aspect.ratio = 1)

# for ea. cluster, sample equal group sizes
cs <- data.table(cd_df, keep.rownames = TRUE) %>%
    split(by = c("cluster_id", "group_id"), 
        flatten = FALSE, sorted = TRUE) %>% 
    map_depth(2, "rn") %>% 
    map_depth(1, function(u) {
        n <- min(sapply(u, length))
        lapply(u, sample, n)
    }) %>% unlist %>% unname
df <- df[sample(cs, length(cs)), ] # randomize rows
for (id in ids) {
    p <- .plot_dr(cd_df, id, get(paste0(id, "_pal")))
    fn <- sprintf(file.path(dir, "magl-umap_%s.pdf"), id)
    ggsave(fn, p, 
        width = 10, height = 10, units = "cm", 
        dpi = 300, useDingbats = FALSE)
}

# pseudobulk-level MDS plot ----------------------------------------------------
# mds <- pbMDS(pb) + scale_color_manual(values = cluster_id_pal)
# xy <- mds$data[, c(1, 2)]
# .get_lims <- function(u, x = 0.25) {
#     x <-  1/x
#     c(floor(min(u)*x)/x, ceiling(max(u)*x)/x)
# }
# lims <- apply(xy, 2, .get_lims)
# p <- mds + .prettify("bw") + theme(legend.position = "none") +
#     scale_shape_manual(values = c("WT" = 17, "LPS" = 8)) +
#     scale_x_continuous(limits = lims[, 1], breaks = seq(-2,2,2), expand = c(0,0.25)) +
#     scale_y_continuous(limits = lims[, 2], breaks = seq(-2,2,2), expand = c(0,0.25))
# p$layers[[1]]$aes_params$size <- 1.5
# ggsave(file.path(dir, "magl-pb_mds.pdf"), p, 
#     width = 5, height = 5, units = "cm", 
#     dpi = 300, useDingbats = FALSE)

# heatmap of cross-group logFCs ------------------------------------------------

# construct SCE containing logFCs for ea. cluster
ref <- pb$group_id == "vehicle"
lfc <- lapply(assays(pb), function(u) {
    y <- DGEList(u)
    y <- calcNormFactors(y)
    lcpm <- log1p(cpm.DGEList(y))
    lfc <- lcpm - rowMeans(lcpm[, ref])
    SingleCellExperiment(lfc, colData = colData(pb))
}) %>% do.call(what = cbind)
lfc$cluster_id <- factor(rep(assayNames(pb), each = ns), levels = kids)

# consensus clustering of genes
cc <- M3C(t(assay(lfc[deg, ])), method = 2)
cc_ids <- cc$realdataresults[[3]]$assignments

cc_cols <- brewer.pal(length(unique(cc_ids)), "Set2")
names(cc_cols) <- unique(cc_ids)
row_anno <- rowAnnotation(
    show_annotation_name = FALSE,
    df = data.frame(consensus_id = cc_ids),
    col = list(consensus_id = cc_cols),
    simple_anno_size = unit(2, "mm"),
    annotation_legend_param = lgd_aes)

qs <- quantile(assay(lfc[deg, ]), c(0.01, 0.99))
brks <- .get_brks(qs, 0.5)
cols <- colorRamp2(brks, hm_cols)

hm <- .plot_hm(lfc[deg, ], cols, brks, 
    col_anno, row_anno, cc_ids, row_title = NULL, 
    cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
    use_raster = TRUE, raster_device = "CairoPNG")
hm@top_annotation@anno_list[[1]]@show_legend <- FALSE
hm@top_annotation@anno_list[[2]]@show_legend <- FALSE
g <- grid.grabExpr(draw(hm))

ggsave(file.path(dir, "magl-lfc_hm.pdf"), grid.draw(g),
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)

# consensus ID 3 - top GO term -------------------------------------------------
go <- lapply(split(names(cc_ids), cc_ids), function(gs) {
    ss1 <- strsplit(gs, ".", fixed = TRUE)
    ss2 <- strsplit(rownames(sce), ".", fixed = TRUE)
    fs1 <- unique(sapply(ss1, .subset, 2))
    fs2 <- sapply(ss2, .subset, 2)
    l <- as.numeric(fs2 %in% fs1)
    names(l) <- fs2
    go <- new("topGOdata", 
        nodeSize = 20,
        ontology = "BP",
        allGenes = l,
        geneSel = function(u) u == 1,
        annot = annFUN.org,
        mapping = "org.Mm.eg.db",
        ID = "symbol")
    res <- runTest(go, algorithm = "classic", statistic = "fisher")
    GenTable(go, classicFisher = res, orderBy = "classicFisher", topNodes = 100)
})

head(go[[3]])

# cell-level linear decomposition ----------------------------------------------
cs_by_k <- split(colnames(sce), sce$cluster_id)
ld <- lapply(kids, function(k) {
    # subset SCE to inlcude
    # - genes that are DE in cluster k
    # - cells that have been assigned to cluster k
    gs <- deg_by_k[[k]]
    cs <- cs_by_k[[k]]
    if (length(gs) == 0 | length(cs) == 0) 
        return(NULL)
    es <- logcounts(sce[gs, cs])
    # decompose values as a linear combination of the group's logFCs
    ld <- t(es) %*% tbl[[k]]$logFC
    # assure values are comparable across clusters
    ld[, 1] / max(ld)
})
rmv <- vapply(ld, is.null, logical(1))
ld <- set_names(unlist(ld[!rmv]), unlist(cs_by_k[!rmv]))
cd_df[names(ld), "ld"] <- ld

# violins by cluster-sample
p <- ggplot(cd_df, 
    aes(x = sample_id, y = ld, fill = group_id)) +
    facet_wrap("cluster_id", nrow = 1) +
    geom_hline(yintercept = 0, size = 0.1, lty = 2) +
    scale_fill_manual(values = gcols) + 
    geom_violin(size = 0.1, width = 1) + 
    geom_boxplot(size = 0.1, width = 0.2, show.legend = FALSE,
        outlier.shape = 16, outlier.size = 0.1) +
    labs(x = NULL, y = "effect coefficient") + 
    scale_y_continuous(limits = c(-0.6, 1.1), breaks = seq(-0.5,1,0.5), expand = c(0,0)) +
    .prettify("bw") + theme(aspect.ratio = 1, 
        legend.position = "none", panel.spacing = unit(1, "mm"), 
        strip.background = element_blank(), strip.text = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggsave(file.path(dir, "magl-ld_violins.pdf"), p,
    width = 13.5, height = 2.5, units = "cm",
    dpi = 300, useDingbats = FALSE)

# split by consensus ID --------------------------------------------------------

cs_by_k <- split(colnames(sce), sce$cluster_id)
ld <- lapply(kids, function(k) {
    gs <- deg_by_k[[k]]
    cs <- cs_by_k[[k]]
    re <- tbl[[k]]
    lapply(unique(cc_ids), function(c) {
        gs <- names(which(cc_ids[gs] == c))
        if (length(gs) == 0 | length(cs) == 0) 
            return(NULL)
        re <- filter(re, gene %in% gs)
        es <- logcounts(sce[gs, cs])
        # decompose values as a linear combination of the group's logFCs
        ld <- t(es) %*% re[, "logFC"]
        # assure values are comparable across clusters
        ld[, 1] / max(ld)
    })
})
lds <- setNames(unlist(ld), unlist(sapply(ld, sapply, names)))
df <- data.frame(ld = lds, cd_df[names(lds), ids],
    cc_id = rep.int(rep(unique(cc_ids), nk), unlist(map_depth(ld, 2, length))))
df <- df[is.finite(df$ld), ]

p <- ggplot(df, aes(x = sample_id, y = ld, fill = group_id)) +
    facet_grid(rows = vars(cc_id), cols = vars(cluster_id), scales = "free_y") +
    geom_hline(yintercept = 0, size = 0.1, lty = 2) +
    scale_fill_manual(values = gcols) + 
    geom_violin(size = 0.1, width = 1) + 
    geom_boxplot(size = 0.1, width = 0.2, show.legend = FALSE,
        outlier.shape = 16, outlier.size = 0.1) +
    labs(x = NULL, y = "effect coefficient") + 
   # scale_y_continuous(limits = c(-0.6, 1.1), breaks = seq(-0.5,1,0.5), expand = c(0,0)) +
    .prettify("bw") + theme(aspect.ratio = 1, 
        legend.position = "bottom", strip.text = element_text(size = 5),
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggsave(file.path(dir, "magl-ld_violins_split.pdf"), p,
    width = 15, height = 6.5, units = "cm",
    dpi = 300, useDingbats = FALSE)

# GSEA -------------------------------------------------------------------------

m_df <- msigdbr(
    species = "Mus musculus") %>%
    dplyr::filter(gs_cat %in% c("H", "C5", "C7"))

ss <- strsplit(rownames(sce), ".", fixed = TRUE)
rowData(sce)$ensembl_id <- sapply(ss, .subset, 1)
rowData(sce)$symbol <- sapply(ss, .subset, 2)

dat <- lapply(res$data, function(u) {
    m <- match(rownames(u), rownames(sce))
    ss <- strsplit(rownames(u), ".", fixed=TRUE)
    u$genes <- as.data.frame(rowData(sce))[m, ]
    return(u)
})

sets <- split(m_df$gene_symbol, m_df$gs_name)
n <- vapply(sets, length, numeric(1))
sets <- sets[n >= 20 & n <= 1000]

mm <- model.matrix(~0+pb$group_id, levels = levels(pb$group_id))
rownames(mm) <- colnames(pb)
colnames(mm) <- levels(pb$group_id)
contrast <- makeContrasts("LPS-vehicle", levels = pb$group_id)

gs_dat <- mapply(function(uu, vv) {
    inds <- ids2indices(sets, uu$genes$symbol, remove.empty = TRUE)
    d <- mm[colnames(uu),]
    v <- voom(uu, d)
    f <- lmFit(v, d)
    f <- eBayes(f)
    cf <- contrasts.fit(f, contrast)
    cf <- eBayes(cf)
    list(indices = inds, voom = v, design = d, 
        cluster_id = vv, contrasts.fit = cf)
}, dat, names(dat), SIMPLIFY = FALSE)

gs_df <- lapply(gs_dat, function(u)
    camera(u$voom, u$indices, u$design, contrast) %>% 
        rownames_to_column("geneset")) %>% 
    bind_rows(.id = "cluster_id")

cats <- gs_df %>% 
    dplyr::filter(FDR < 1e-20) %>%
    pull(geneset) %>% unique
length(cats)
mat <- gs_df %>% 
    dplyr::filter(geneset %in% cats) %>%
    dplyr::mutate(neg_log10_fdr = -log10(FDR)) %>% 
    reshape2::acast(cluster_id ~ geneset, value.var = "neg_log10_fdr") %>% 
    set_colnames(gsub("/*([^_]*)_(.*)", "\\2", colnames(.))) %>% 
    set_colnames(strtrim(colnames(.), 30)) 

hm <- Heatmap(mat, 
    #name = ,
    col = circlize::colorRamp2(c(0, 10, 20, 40), 
    c("white", "cornflowerblue", "violet", "red")),  
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 4),
    heatmap_legend_param = list(
        direction = "horizontal",
        title_position = "lefttop",
        title = expression("-log"[10](FDR)),
        labels_gp = gpar(fontsize = 6)),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"))

ggsave(file.path(dir, "magl-gsea_hm.pdf"), 
    grid.draw(grid.grabExpr(draw(hm, heatmap_legend_side = "bottom"))),
    width = 15, height = 10, units = "cm",
    dpi = 300, useDingbats = FALSE)

# consensus ID 3 ---------------------------------------------------------------
gs <- names(which(cc_ids == 3))
gs <- intersect(gs, unique(unlist(.deg(res, 1e-10))))

qs <- quantile(assay(lfc[gs, ]), c(0.01, 0.99))
brks <- .get_brks(qs, 0.5)
cols <- colorRamp2(brks, hm_cols)

hm <- .plot_hm(lfc[gs, ], cols, brks, 
    col_anno, row_title = NULL, row_nms = TRUE,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
    row_names_gp = gpar(fontsize = 5),
    column_title_gp = gpar(fontsize = 8),
    use_raster = TRUE, raster_device = "CairoPNG")

ggsave(file.path(dir, "magl-lfc_cc3.pdf"), 
    grid.draw(grid.grabExpr(draw(hm))),
    width = 15, height = 12, units = "cm",
    dpi = 300, useDingbats = FALSE) 

# GO_CHEMOKINE_RECEPTOR_BINDING ------------------------------------------------

fs <- sets[["GO_CHEMOKINE_RECEPTOR_BINDING"]]
idx <- match(fs, rowData(sce)$symbol, nomatch = 0)
m <- match(rownames(sce)[idx], rownames(lfc))

qs <- quantile(assay(lfc[m, ]), c(0.01, 0.99))
brks <- .get_brks(qs, 0.5)
cols <- colorRamp2(brks, hm_cols)

hm <- .plot_hm(lfc[m, ], cols, brks, 
    col_anno, row_title = NULL, row_nms = TRUE,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
    row_names_gp = gpar(fontsize = 5),
    column_title_gp = gpar(fontsize = 8),
    use_raster = TRUE, raster_device = "CairoPNG")

ggsave(file.path(dir, "magl-lfc_chemokine.pdf"), 
    grid.draw(grid.grabExpr(draw(hm))),
    width = 15, height = 7, units = "cm",
    dpi = 300, useDingbats = FALSE)


# interferon_gama_response XOR interferon_alpha_response -----------------------
ifs <- sapply(c("interferon_alpha", "interferon_gamma"), 
    grep, names(sets), ignore.case = TRUE) %>% 
    sapply(function(u) unlist(sets[u]))
fs <- setdiff(unlist(ifs), intersect(ifs[[1]], ifs[[2]]))

m0 <- match(fs, rowData(sce[deg, ])$symbol, nomatch = 0)
m <- match(rownames(lfc), deg[m0], nomatch = 0)

qs <- quantile(assay(lfc[m, ]), c(0.01, 0.99))
brks <- .get_brks(qs, 1)
cols <- colorRamp2(brks, hm_cols)

w <- as.numeric(fs[m0 != 0] %in% ifs[[1]]) + 1
row_cols <- setNames(c("skyblue", "lightgreen"), c("IFN-alpha", "IFN-beta"))
row_anno <- rowAnnotation(
    show_annotation_name = FALSE,
    df = data.frame(response = names(row_cols)[w]),
    col = list(response = row_cols),
    simple_anno_size = unit(1, "mm"),
    annotation_legend_param = c(lgd_aes))

hm <- .plot_hm(lfc[m, ], cols, brks, 
    col_anno, row_anno, row_split = w, row_title = NULL, row_nms = TRUE,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
    row_names_gp = gpar(fontsize = 5),
    column_title_gp = gpar(fontsize = 8),
    use_raster = TRUE, raster_device = "CairoPNG")

ggsave(file.path(dir, "magl-lfc_inf.pdf"), 
    grid.draw(grid.grabExpr(draw(hm))),
    width = 15, height = 12, units = "cm",
    dpi = 300, useDingbats = FALSE)

