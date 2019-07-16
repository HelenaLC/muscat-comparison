# load packages
suppressMessages({
    library(muscat)
    library(SingleCellExperiment)
})

# load data
sce <- readRDS(file.path(snakemake@config$raw_data, "ref_magl.rds"))

# keep reference samples only
sce <- sce[, sce$group_id == "WT"]

# prep. SCE for 'muscat'
sce <- prepSCE(sce, 
    cluster_id = "cluster_id",
    sample_id = "sample_id",
    group_id = "group_Id",
    drop = TRUE)

# remove slots other than counts
assays(sce) <- SimpleList(counts = counts(sce))

# remove dimension reductions
reducedDims(sce) <- NULL

# write SCE to .rds
saveRDS(sce, "/Users/helena/Desktop/sce0_magl.rds")
