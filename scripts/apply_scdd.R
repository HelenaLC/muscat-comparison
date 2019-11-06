suppressMessages({
    library(dplyr)
    library(scater)
    library(scDD)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_scdd <- function(sce, pars, ds_only = TRUE) {
    t <- system.time({
        if (ds_only) {
            normcounts(sce) <- switch(pars$assay,
                logcounts = 2^logcounts(sce)-1,
                vstresiduals = exp(assay(sce, "vstresiduals")))
        } else {
            normcounts(sce) <- switch(pars$assay, 
                logcounts = normcounts(logNormCounts(computeLibraryFactors(sce), log = FALSE)),
                vstresiduals = exp(vst(counts(sce), show_progress = FALSE)$y))
        }
        res <- tryCatch(run_scdd(sce), error = function(e) e)
    })[[3]]
    list(rt = t, tbl = res)
}

run_scdd <- function(sce) {
    kids <- levels(sce$cluster_id)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    if (is(normcounts(sce), "dgCMatrix"))
        normcounts(sce) <- as.matrix(normcounts(sce))
    suppressMessages(
        res <- lapply(kids, function(k) {
            res <- results(scDD(sce[, cells_by_k[[k]]], 
                min.nonzero = 20, condition = "group_id",
                categorize = FALSE, testZeroes = FALSE,
                param = BiocParallel::MulticoreParam(workers = 1)))
            data.frame(
                gene = rownames(sce), 
                cluster_id = k,
                p_val = res$nonzero.pvalue, 
                p_adj.loc = res$nonzero.pvalue.adj,
                row.names = NULL, 
                stringsAsFactors = FALSE)
        }) 
    )
    df <- bind_rows(res)
    df$p_adj.glb <- p.adjust(df$p_val)
    return(df)
}
