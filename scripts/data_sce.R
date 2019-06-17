# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(SingleCellExperiment)
    }))

# get data parameters
md <- readRDS(snakemake@input$data_pars)

# load data
if (length(grep(".rds", md$dir_raw)) == 1) {
    x <- readRDS(md$dir_raw)
    rd <- rowData(x)
    cd <- colData(x)
    x <- counts(x)
} else {
    x <- lapply(md$dir_raw, Matrix::readMM)
    x <- do.call("cbind", x)
    x <- as(x, "dgCMatrix") 
}

# load metadata
if (length(md$dir_cell_md) != 0) {
    cd <- read.delim(md$dir_cell_md, stringsAsFactors = FALSE)
    colnames(x) <- rownames(cd) <- paste0("cell", seq_len(nrow(cd)))
}
if (length(md$dir_gene_md) != 0) {
    rd <- read.delim(md$dir_gene_md, stringsAsFactors = FALSE, 
        header = FALSE, col.names = c("gene", "feature"))
    rownames(x) <- rownames(rd) <- rd$gene
}

# keep cells with at least 200 expressed genes,
# and genes detected in at least 100 cells
cs_keep <- Matrix::colSums(x > 0) >= 200
gs_keep <- Matrix::rowSums(x > 0) >= 100

x <- x[gs_keep, cs_keep]
cd <- cd[colnames(x), , drop = FALSE]
rd <- rd[rownames(x), , drop = FALSE]

if (!is.null(md$filter_cells)) {
    cells_keep <- md$filter_cells(cd)
    cd <- cd[cells_keep, ]
    x <- x[, rownames(cd)]
}

# construct colData
for (i in c("sample_id", "group_id", "cluster_id")) {
    if (is.function(md[[i]])) {
        ids <- md[[i]](cd)
    } else {
        ids <- cd[[md[[i]]]]
    }
    assign(i, ids)
}
cd <- data.frame(sample_id, group_id, cluster_id)

# construct SingleCellExperiment
sce <- SingleCellExperiment(
    assays = list(counts = x),
    rowData = rd, colData  = cd)

ei <- muscat:::.make_ei(sce)
metadata(sce)$experiment_info <- ei

# write to .rds
saveRDS(sce, snakemake@output[[1]])