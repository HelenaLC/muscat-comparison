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
            error = function(e) NULL, 
            run_scdd(sce))
    })[[3]]
    
    # return results
    list(rt = t, tbl = res)
}

run_scdd <- function(sce) {
    priors <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
    sce$group_id <- as.numeric(sce$group_id)
    kids <- levels(sce$cluster_id)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    suppressMessages(
        lapply(kids, function(k) {
            res <- results(scDD(sce[, cells_by_k[[k]]],
                condition = "group_id", prior_param = priors, 
                testZeroes = FALSE, min.nonzero = 20))
            data.frame(
                gene = rownames(sce),
                cluster_id = k,
                p_val = res$nonzero.pvalue, 
                p_adj.loc = res$nonzero.pvalue.adj,
                row.names = NULL,
                stringsAsFactors = FALSE)
        }) %>% bind_rows %>% dplyr::mutate(p_adj.glb = p.adjust(p_val))
    )
}
