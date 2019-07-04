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
                    normcounts = normcounts(normalize(sce, return_log = FALSE)),
                    logcounts = logcounts(normalize(sce)),
                    vstcounts = vst(counts(sce), show_progress = FALSE)$y,
                    cpm = calculateCPM(counts(sce)),
                    logcpm = log2(calculateCPM(counts(sce)) + 1))))
        }
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]

    # run & time DS analysis
    t2 <- system.time(suppressWarnings(
        res <- tryCatch(error = function(e) NULL, 
            pbDS(pb, method = pars$method, verbose = FALSE))))[[3]]

    # return results
    list(rt = c(t1, t2), tbl = dplyr::bind_rows(res$table[[1]]))
}