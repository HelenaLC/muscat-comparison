suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    # run & time aggregation
    t1 <- system.time({
        a <- pars$assay
        if (!ds_only) {
            suppressWarnings(suppressMessages(
                assays(sce)[[a]] <- switch(a, 
                    counts = counts(sce),
                    cpm = calculateCPM(counts(sce)),
                    logcounts = logNormCounts(computeLibraryFactors(sce)),
                    vstresiduals = vst(counts(sce), show_progress = FALSE)$y)))
        }
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]

    # run & time DS analysis
    t2 <- system.time({
        res <- tryCatch(
            pbDS(pb, method = pars$method, verbose = FALSE),
            error = function(e) e)
        if (!inherits(res, "error"))
            res <- dplyr::bind_rows(res$table[[1]])
    })[[3]]

    # return results
    list(rt = c(t1, t2), tbl = res)
}
