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

# load data parameters
md <- readRDS(snakemake@input$data_pars)

# load data
sce <- readRDS(snakemake@input$dir_data)

# filter cluster-samples with >= 100 cells
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, n = 100)
sce <- .filter_sce(sce, dimnames(n_cells))

# remove group 2 samples
cs_keep <- sce$group_id == md$group_keep
sce <- sce[, cs_keep]

# update colData & metdata
cd <- data.frame(colData(sce))
cd <- mutate_if(cd, is.factor, droplevels)
colData(sce) <- S4Vectors::DataFrame(cd, row.names = colnames(sce))
metadata(sce)$experiment_info <- muscat:::.make_ei(sce)

# fit NB
y <- DGEList(counts(sce))
mm <- model.matrix(~ 0 + sce$sample_id)
y <- estimateDisp(y, mm)
y <- glmFit(y)

# update row- & colData
sce$offset <- c(y$offset)
rowData(sce)$dispersion <- y$dispersion
rowData(sce) <- DataFrame(rowData(sce), y$coefficients %>% 
        set_colnames(paste("beta", levels(sce$sample_id), sep = ".")))

# write to .rds
saveRDS(sce, snakemake@output[[1]])