source(".Rprofile")
suppressMessages({
    library(dplyr)
    library(edgeR)
    library(Matrix)
    library(SingleCellExperiment)  
})
#-------------------------------------------------------------------------------
sce <- readRDS(snakemake@input$sce)

# keep genes w/ count > 1 in at least 10 cells;
# keep cells w/ at least 100 detected genes
gs <- rowSums(counts(sce) > 1) >= 10
cs <- colSums(counts(sce) > 0) >= 100
sce <- sce[gs, cs]

# keep cluster-samples w/ at least 100 cells
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, n = 100)

kids <- rownames(n_cells)
sids <- colnames(n_cells)
sce <- .filter_sce(sce, kids, sids)

# esimate/compute gene/cell parameters
y <- DGEList(counts(sce))
mm <- model.matrix(~ 0 + sce$sample_id)
y <- estimateDisp(y, mm)
y <- glmFit(y, prior.count = 0)

# update row- & colData
sce$offset <- c(y$offset)
rowData(sce)$dispersion <- y$dispersion
betas <- paste("beta", levels(sce$sample_id), sep = ".")
colnames(y$coefficients) <- betas
rowData(sce)[, betas] <- y$coefficients

# write SCE to .rds
saveRDS(sce, snakemake@output$sce)