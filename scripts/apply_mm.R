apply_mm <- function(x, pars) {
    # get method parameters
    args <- as.list(args("mmDS"))
    pars <- pars[names(pars) %in% names(args)]
    
    # run & time method
    t <- system.time({
        if (pars$covs == "dr") {
            colData(x)$dr <- Matrix::colMeans(counts(x) > 0)
        } else {
            pars$covs <- NULL
        }
        res <- tryCatch(error = function(e) NULL, {
            do.call(mmDS, c(list(x, verbose = FALSE), pars))
        })
    })[[3]]

    # return results
    list(rt = t, res = bind_rows(res))
}
