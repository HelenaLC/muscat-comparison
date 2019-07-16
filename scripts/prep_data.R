# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(edgeR)
        library(magrittr)
        library(SingleCellExperiment)
    }))

# source utils
source(snakemake@config$utils)

# load data
sce <- readRDS(snakemake@input$sce)

# filter cluster-samples with >= 100 cells
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, n = 100)

kids <- rownames(n_cells)
sids <- colnames(n_cells)
sce <- .filter_sce(sce, kids, sids)

# remove group 2 samples
cs_keep <- sce$group_id == md$group_keep
sce <- sce[, cs_keep]

# update colData & metdata
cd <- data.frame(colData(sce))
cd <- mutate_if(cd, is.factor, droplevels)
colData(sce) <- DataFrame(cd, row.names = colnames(sce))
metadata(sce)$experiment_info <- muscat:::.make_ei(sce)

# fit NB
y <- DGEList(counts(sce))
mm <- model.matrix(~ 0 + sce$sample_id)
y <- estimateDisp(y, mm)
y <- glmFit(y, prior.count = 0)

# update row- & colData
sce$offset <- c(y$offset)
rowData(sce)$dispersion <- y$dispersion

betas <- paste("beta", levels(sce$sample_id), sep = ".")
coefs <- set_colnames(y$coefficients, betas)
rowData(sce) <- DataFrame(rowData(sce), coefs)

# write to .rds
saveRDS(sce, snakemake@output[[1]])