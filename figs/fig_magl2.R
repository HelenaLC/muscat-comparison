# marker genes heatmap /dotplot ------------------------------------------------
known_markers <- list(
    astrocytes = c("Aqp4", "Gfap", "Fgfr3"),
    endothelial = c("Cldn5","Nostrin"),
    microglia = c("C1qb","Tyrobp"),
    neuronal = c("Snap25", "Stmn2"),
    neuronal_excitatory = "Slc17a7",
    neuronal_inhibitory = "Gad1",
    oligodendrocytes = "Opalin",
    OPC = "Pdgfra")

# grep full names
known_markers <- lapply(known_markers, sapply, function(g) 
    grep(paste0(g, "$"), rownames(sce), value = TRUE))

# make labels
gs <- gsub(".*\\.", "", unlist(known_markers))
ks <- rep.int(names(known_markers), vapply(known_markers, length, numeric(1)))
labs <- sprintf("%s(%s)", gs, ks)

scran_markers <- findMarkers(sce, 
    clusters = sce$cluster_id, block = sce$sample_id, 
    direction = "up", lfc = 2, full.stats = TRUE)

# pull top
top <- lapply(scran_markers, function(u)
    data.frame(gene = rownames(u), u, check.names = FALSE)) %>% 
    bind_rows(.id = "cluster_id") %>% 
    filter(Top == 1)
scran_gs <- as.character(top$gene)

# split cells by cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)

# compute cluster-marker means
markers <- c(scran_gs, unlist(known_markers))
ms_by_k <- vapply(cs_by_k, function(i)
    Matrix::rowMeans(logcounts(sce)[markers, i, drop = FALSE]),
    numeric(length(markers)))

mat <- muscat:::.scale(ms_by_k)
fs <- gsub("^.*\\.", "", rownames(sce))
m <- match(markers, rownames(sce))
rownames(mat) <- fs[m]

u <- names(which(table(rownames(mat)) == 1))

cols <- setNames(
    muscat:::.cluster_colors[seq_len(nlevels(sce$cluster_id) + length(known_markers))], 
    c(levels(sce$cluster_id), names(known_markers)))
row_anno <- rowAnnotation(
    df = data.frame(label = c(top$cluster_id, ks)),
    col = list(label = cols),
    gp = gpar(col = "white")) 

Heatmap(mat[u, ],
    name = "scaled avg.\nexpression",
    col = viridis(10),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = "cluster_id",
    column_title_side = "bottom",
    rect_gp = gpar(col = "white"),
    #left_annotation = row_anno,
    row_names_side = "left")


means_by_cluster <- lapply(known_markers, function(gs)
    vapply(cs_by_k, function(i)
        Matrix::rowMeans(logcounts(sce)[gs, i, drop = FALSE]), 
        numeric(length(gs))))

# prep. for plotting & scale b/w 0 and 1
mat <- do.call("rbind", means_by_cluster)
mat <- muscat:::.scale(mat)
rownames(mat) <- gs

cols <- muscat:::.cluster_colors[seq_along(known_markers)]
cols <- setNames(cols, names(known_markers))

Heatmap(mat,
    name = "scaled avg.\nexpression",
    col = viridis(10),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_title = "cluster_id",
    column_title_side = "bottom",
    rect_gp = gpar(col = "white"),
    left_annotation = row_anno)