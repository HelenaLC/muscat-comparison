set.seed(2903)
print(args)
suppressMessages({
    library(countsimQC)
    library(DESeq2)
    library(muscat)
    library(SummarizedExperiment)
})

sce <- readRDS(args$sce)

sim <- simData(sce,
    nk = 3, ns = 3,
    ng = nrow(sce), 
    nc = ncol(sce), 
    p_dd = diag(6)[1, ], 
    probs = list(NULL, NULL, c(1, 0)))

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
    outputFile = basename(args$html),
    outputDir = dirname(args$html),
    outputFormat = "html_document",
    maxNForCorr = 200,
    maxNForDisp = 500,
    forceOverwrite = TRUE)
