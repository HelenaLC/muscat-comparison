suppressMessages({
    library(muscat)
    library(SingleCellExperiment)
})

apply_mm <- function(sce, pars, ds_only = TRUE) {
    pars <- pars[names(pars) != "id"]
    if (is.null(pars$n_threads))
        pars$n_threads <- 1
    t <- system.time({
        pars[sapply(pars, `==`, "")] <- NULL
        res <- tryCatch(
            do.call(mmDS, c(list(sce, verbose = FALSE), pars)),
            error = function(e) e)
        if (!inherits(res, "error"))
            res <- dplyr::bind_rows(res)
    })[[3]]
    list(rt = t, tbl = res)
}
