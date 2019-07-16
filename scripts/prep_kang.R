# load data
sce <- readRDS(snakemake@input$sce)

# remove multiplets & unassignet cells
sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]

# keep reference samples only
sce <- sce[, sce$stim == "ctrl"]

# construct sample IDs
sce$sample_id <- factor(paste0(sce$stim, sce$ind))

# prep. SCE for 'muscat'
sce <- muscat::prepSCE(sce, 
    cluster_id = "cell",
    sample_id = "sample_id",
    group_id = "stim",
    drop = TRUE)

# remove slots other than counts
assays(sce) <- S4Vectors::SimpleList(counts = counts(sce))

# remove dimension reductions
reducedDims(sce) <- NULL

# write SCE to .rds
saveRDS(sce, "/Users/helena/Desktop/sce0_kang.rds")
