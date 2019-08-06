options(conflicts.policy = list(warn = FALSE))
suppressMessages(library(SingleCellExperiment))

# load data & method parameters
sce <- readRDS(snakemake@input$sce)
sce$cluster_id <- factor(sce$cluster_id)
meth_pars <- as.list(jsonlite::fromJSON(snakemake@input$meth_pars))

# get method wrapper 'apply_xx()'
source(fun <- snakemake@input$fun)
fun <- gsub("(.R)", "", basename(fun))

if (isTRUE(grep("MM", meth_pars$id) == 1))
    meth_pars$n_threads <- 10

# run method & write results to .rds
res <- get(fun)(sce, meth_pars, ds_only = FALSE)
saveRDS(res, snakemake@output$res)
