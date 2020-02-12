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

