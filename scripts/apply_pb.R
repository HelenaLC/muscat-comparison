suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    t1 <- system.time({
        a <- pars$assay
        if (!ds_only) 
            assay(sce, a) <- switch(a, 
                counts = counts(sce),
                cpm = calculateCPM(counts(sce)),
                logcounts = logNormCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]
    t2 <- system.time({
        res <- tryCatch(
            pbDS(pb, method = pars$method, filter = "genes", verbose = FALSE),
            error = function(e) e)
        if (!inherits(res, "error"))
            res <- dplyr::bind_rows(res$table[[1]])
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}