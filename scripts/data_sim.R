# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(edgeR)
        library(magrittr)
        library(S4Vectors)
        library(SingleCellExperiment)
    }))

# source utils
source(snakemake@config$utils)

# load data parameters
md <- readRDS(snakemake@input$data_pars)

# load data
sce <- readRDS(snakemake@input$dir_data)

# remove group 2 samples
cs_keep <- sce$group_id == md$group_keep
sce <- sce[, cs_keep]
sce <- .update_cd(sce)

# filter cluster-samples with >= 100 cells
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, n = 100)
sce <- .filter_sce(sce, dimnames(n_cells))
sce <- .update_cd(sce)

# update metdata
metadata(sce)$experiment_info <- muscat:::.make_ei(sce)

# fit NB
y <- DGEList(counts(sce))
y <- calcNormFactors(y)
mm <- model.matrix(~ 0 + sce$sample_id)
y <- estimateDisp(y, mm)
y <- glmFit(y, prior.count = 0)

# update row- & colData
sce$offset <- c(y$offset)
rowData(sce)$dispersion <- y$dispersion
bs <- paste("beta", levels(sce$sample_id), sep = ".")
rowData(sce) <- DataFrame(rowData(sce), set_colnames(y$coefficients, bs))

# write to .rds
saveRDS(sce, snakemake@output[[1]])