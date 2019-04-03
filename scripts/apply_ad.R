apply_ad <- function(x, pars) {
    # get method parameters
    args <- as.list(args("run_ad"))
    pars <- pars[names(pars) %in% names(args)]
    
    # run & time method
    t <- system.time(
        res <- tryCatch(error = function(e) NULL, 
            do.call(run_ad, c(list(x), pars))))[[3]]
    
    # return results
    list(rt = t, res = res)
}


# ------------------------------------------------------------------------------
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(kSamples)
        library(purrr)
        library(SingleCellExperiment)
    })
)

run_ad <- function(x) {
    nk <- length(kids <- levels(x$cluster_id))
    p_val <- apply(logcounts(x), 1, function(y) {
        re <- ad.test.combined(
            y ~ group_id | cluster_id, 
            data = data.frame(y, colData(x)))
        map(re$ad.list, 5)
    }) %>% map(unlist) %>% unlist
    data.frame(
        stringsAsFactors = FALSE,
        gene = rep(rownames(x), each = nk),
        cluster_id = rep(kids, nrow(x)),
        p_val, p_adj.glb = p.adjust(p_val)) %>% 
        group_by(cluster_id) %>% 
        mutate(p_adj.loc = p.adjust(p_val)) %>% 
        data.frame
}

