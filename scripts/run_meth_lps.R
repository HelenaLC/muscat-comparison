suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(SingleCellExperiment)
})

# wcs <- list(mid = "edgeR.sum.counts")
# args <- list(
#     fun = "scripts/apply_pb.R",
#     sce = "LPS/output/SCE_annotation.rds",
#     meth_pars = "meta/meth_pars/edgeR.sum.counts.json")

sce <- readRDS(args$sce)
nk <- length(kids <- levels(sce$cluster_id))
meth_pars <- as.list(fromJSON(args$meth_pars))

# increase number of threads used by mixed model methods
if (grepl("MM-", wcs$mid)) meth_pars$n_threads <- 30

# run method & write results to .rds
source(fun <- args$fun)
fun <- gsub("(.R)", "", basename(fun))
res <- get(fun)(sce, meth_pars, ds_only = FALSE)$tbl

# assure all gene-cluster combinations are presents in results table
if (!inherits(res, "error")) {
    res <- mutate_at(res, c("gene", "cluster_id"), as.character)
    res <- as.data.frame(rowData(sce)) %>% 
        mutate(gene = rownames(.)) %>% 
        replicate(n = nk, simplify = FALSE) %>% 
        bind_rows %>% 
        mutate(cluster_id = rep(kids, each = nrow(sce))) %>% 
        left_join(res, by = c("gene", "cluster_id")) %>% 
        mutate(mid = wcs$mid)
}

saveRDS(res, args$res)