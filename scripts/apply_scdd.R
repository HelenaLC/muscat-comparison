apply_scdd <- function(sce, pars) {
    # get method parameters
    args <- as.list(args("run_scdd"))
    pars <- pars[names(pars) %in% names(args)]
    
    # run & time method
    t <- system.time(
        res <- tryCatch(error = function(e) NULL, 
            do.call(run_scdd, c(list(sce), pars))))[[3]]
    
    # return results
    list(rt = t, res = res)
}

# ------------------------------------------------------------------------------
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(scDD)
    })
)

run_scdd <- function(x) {
    priors <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
    x$group_id <- as.numeric(x$group_id)
    kids <- levels(x$cluster_id)
    cells_by_k <- split(colnames(x), x$cluster_id)
    suppressMessages(
        lapply(kids, function(k) {
            res <- results(scDD(x[, cells_by_k[[k]]],
                condition = "group_id", prior_param = priors, 
                testZeroes = FALSE, min.size = 3, min.nonzero = NULL,
                param = BiocParallel::MulticoreParam(workers = 1)))
            data.frame(
                gene = rownames(x),
                cluster_id = k,
                p_val = res$nonzero.pvalue, 
                p_adj.loc = res$nonzero.pvalue.adj,
                row.names = NULL,
                stringsAsFactors = FALSE)
        }) %>% bind_rows %>% mutate(p_adj.glb = p.adjust(p_val))
    )
}