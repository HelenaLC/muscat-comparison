suppressMessages({
    library(BiocParallel)
    library(dplyr)
    library(MAST)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_mast <- function(sce, pars, ds_only = TRUE) {
    t <- system.time({
        if (ds_only) {
            a <- assay(sce, pars$assay)
        } else {
            a <- switch(pars$assay, 
                logcounts = normalizeCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sim), show_progress = FALSE)$y)
        }
        if (!is.matrix(a)) 
            a <- as.matrix(a)
        assay(sce, "lCount") <- a
        if (is.null(n <- pars$n_threads)) n <- 1
        res <- tryCatch(run_mast(sce, n), error = function(e) e)
    })[[3]]
    list(rt = t, tbl = res)
}

run_mast <- function(sce, n_threads) {
    colData(sce)$wellKey <- colnames(sce)
    rowData(sce)$primerid <- rownames(sce)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    kids <- levels(sce$cluster_id)
    bplapply(kids, function(k) {
        y <- sce[, cells_by_k[[k]]]
        sca <- SceToSingleCellAssay(y, check_sanity = FALSE)
        fit <- suppressMessages(zlm(~ group_id, sca))
        lrt <- suppressMessages(lrTest(fit, "group_id"))
        p_val <- lrt[, "hurdle", "Pr(>Chisq)"]
        data.frame(
            gene = rownames(y), cluster_id = k, 
            p_val, 
            p_adj.loc = p.adjust(p_val), 
            row.names = NULL, 
            stringsAsFactors = FALSE)
    }, BPPARAM = MulticoreParam(n_threads)) %>% 
        bind_rows %>% mutate(p_adj.glb = p.adjust(p_val))
}