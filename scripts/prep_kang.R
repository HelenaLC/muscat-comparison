source(".Rprofile")
suppressMessages({
    library(muscat)
    library(SingleCellExperiment)
})
#-------------------------------------------------------------------------------
sce <- readRDS(snakemake@input$sce) # load data
sce <- sce[, sce$multiplets == "singlet"] # remove multiplets
sce <- sce[, !is.na(sce$cell)] # remove unassigned cells
sce <- sce[, sce$stim == "ctrl"] # keep control samples only
sce$sample_id <- factor(paste0(sce$stim, sce$ind)) # construct sample IDs
sce <- prepSCE(sce, "cell", "sample_id", "stim", TRUE) # prep. SCE for `muscat`
reducedDims(sce) <- NULL # remove dimensionality reductions
saveRDS(sce, snakemake@output$sce) # write SCE to .rds
