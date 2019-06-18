suppressMessages({
    library(dplyr)
    library(scater)
    library(scDD)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_scdd <- function(sce, pars, ds_only = TRUE) {
    # run & time method
    t <- system.time({
        if (!ds_only) {
            suppressWarnings(suppressMessages(
                assays(sce)$normcounts <- switch(pars$assay, 
                    logcounts = normcounts(normalize(sce, return_log = FALSE)),
                    vstcounts = exp(vst(counts(sce), show_progress = FALSE)$y))))
        } else {
            assays(sce)$normcounts <- switch(pars$assay,
                logcounts = 2^normcounts(sce)-1,
                vstcounts = exp(assays(sce)$vstcounts))
        }
        res <- tryCatch(
            error = function(e) e, 
            run_scdd(sce))
    })[[3]]
    
    # return results
    list(rt = t, tbl = res)
}

run_scdd <- function(sce) {
    kids <- levels(sce$cluster_id)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    suppressMessages(
        lapply(kids, function(k) {
            res <- results(scDD(sce[, cells_by_k[[k]]], 
                min.nonzero = 20, condition = "group_id",
                categorize = FALSE, testZeroes = FALSE,
                param = BiocParallel::MulticoreParam(workers = 1)))
            data.frame(
                row.names = NULL, 
                stringsAsFactors = FALSE,
                gene = rownames(sce), 
                cluster_id = k,
                p_val = res$nonzero.pvalue, 
                p_adj.loc = res$nonzero.pvalue.adj)
        }) %>% bind_rows %>% dplyr::mutate(p_adj.glb = p.adjust(p_val))
    )
}
