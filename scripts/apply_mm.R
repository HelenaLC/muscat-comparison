suppressMessages({
    library(muscat)
    library(SingleCellExperiment)
})

apply_mm <- function(sce, pars, ds_only = TRUE) {
    pars <- pars[names(pars) != "id"]
    
    # run & time method
    t <- system.time({
        if (!ds_only & pars$covs == "dr")
            sce$dr <- colMeans(counts(sce) > 0)
        pars[sapply(pars, `==`, "")] <- NULL
        res <- tryCatch(
            do.call(mmDS, c(list(sce, n_threads = 1, verbose = FALSE), pars)),
            error = function(e) e)
        if (!inherits(res, "error"))
            res <- dplyr::bind_rows(res)
    })[[3]]
    
    # return results
    list(rt = t, tbl = res)
}
