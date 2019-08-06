options(conflicts.policy = list(warn = FALSE))
suppressMessages(library(SingleCellExperiment))

# load data & method parameters
sce <- readRDS(snakemake@input$sce)
meth_pars <- as.list(jsonlite::fromJSON(snakemake@input$meth_pars))

# compute new assay slot if required
if (meth_pars$id == "edgeR.sum(scalecpm)")
    assays(sce)$cpm <- scater::calculateCPM(sce)

if (isTRUE(grep("vstresiduals", meth_pars$id) == 1))
    assays(sce)$vstresiduals <- suppressWarnings(
        sctransform::vst(counts(sce), show_progress = FALSE)$y)

# get method wrapper 'apply_xx()'
source(fun <- snakemake@input$fun)
fun <- gsub("(.R)", "", basename(fun))

# run method & write results to .rds
res <- get(fun)(sce, meth_pars)
saveRDS(res, snakemake@output$res)
