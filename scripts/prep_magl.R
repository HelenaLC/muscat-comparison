args <- R.utils::commandArgs(
    trailingOnly = TRUE, 
    asValues = TRUE)

suppressMessages({
    library(muscat)
    library(S4Vectors)
    library(SingleCellExperiment)
})

sce <- readRDS(args$input_sce) # load data
reducedDims(sce) <- NULL # remove dimension reductions
assays(sce) <- SimpleList(counts = counts(sce)) # remove slots other than counts
sce <- sce[, sce$group_id == "WT"] # keep reference samples only
sce <- prepSCE(sce, "cluster_id", "sample_id", "group_id", TRUE) # prep. SCE for 'muscat'
saveRDS(sce, args$output_sce) # write SCE to .rds
