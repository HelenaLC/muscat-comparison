# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(countsimQC)
        library(DESeq2)
        library(muscat)
        library(SummarizedExperiment)
    }))

# load data
sce <- readRDS(snakemake@input$data)
probs <- list(
    table(sce$cluster_id) / ncol(sce),
    table(sce$sample_id)  / ncol(sce),
    NULL)

# simulate data
set.seed(1)
sim <- simData(sce,
    n_genes = nrow(sce),
    n_cells = ncol(sce), 
    p_dd = diag(6)[1, ],
    probs = probs)

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
    outputFile = snakemake@output$res,
    outputDir = ".",
    outputFormat = "html_document",
    maxNForCorr = 1e3,
    maxNForDisp = 2e3,
    forceOverwrite = TRUE)