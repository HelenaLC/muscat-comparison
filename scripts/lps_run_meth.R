options(conflicts.policy = list(warn = FALSE))
suppressMessages(library(SingleCellExperiment))

# load data & method parameters
sce <- readRDS(snakemake@input$sce)
meth_pars <- as.list(jsonlite::fromJSON(snakemake@input$meth_pars))

# get method wrapper 'apply_xx()'
source(fun <- snakemake@input$fun)
fun <- gsub("(.R)", "", basename(fun))

# run method & write results to .rds
res <- get(fun)(sce, meth_pars, ds_only = FALSE)
if (!inherits(res$tbl, "error")) {
    foo <- data.frame(
        gene = rep(rownames(sce), nk),
        cluster_id = rep(kids, each = nrow(sce)),
        method_id = meth_pars$id,
        stringsAsFactors = FALSE)
    res$tbl <- left_join(foo, res$tbl, 
        by = c("gene", "cluster_id"))
}
saveRDS(res, snakemake@output$res)
