source(".Rprofile")
suppressMessages({
    library(muscat)
    library(S4Vectors)
    library(SingleCellExperiment)
})
#-------------------------------------------------------------------------------
sce <- readRDS(snakemake@input$sce) # load data
sce <- sce[, sce$group_id == "WT"] # keep reference samples only
assays(sce) <- SimpleList(counts = counts(sce)) # remove slots other than counts
sce <- prepSCE(sce, "cluster_id", "sample_id", "group_id", TRUE) # prep. SCE for 'muscat'
reducedDims(sce) <- NULL # remove dimension reductions
saveRDS(sce, snakemake@output$sce) # write SCE to .rds
