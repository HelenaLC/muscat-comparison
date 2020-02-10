suppressMessages({
    library(countsimQC)
    library(DESeq2)
    library(muscat)
    library(SummarizedExperiment)
})

sce <- readRDS(sce_fn)
probs <- list(
    table(sce$cluster_id) / ncol(sce),
    table(sce$sample_id)  / ncol(sce),
    NULL)

# simulate data
set.seed(2903)
sim <- simData(args$sce,
    ng = nrow(sce), nc = ncol(sce), 
    p_dd = diag(6)[1, ], probs = probs)

dds_sim <- DESeqDataSetFromMatrix(
    countData = counts(sim),
    colData = colData(sim),
    design = ~ sample_id)

dds_sce <- DESeqDataSetFromMatrix(
    countData = counts(sce),
    colData = colData(sce),
    design = ~ sample_id)

dds_list <- list(
    sim = dds_sim,
    data = dds_sce)

countsimQCReport(
    ddsList = dds_list,
    outputFile = basename(args$report),
    outputDir = dirname(args$report),
    outputFormat = "html_document",
    maxNForCorr = 200,
    maxNForDisp = 500,
    forceOverwrite = TRUE)
