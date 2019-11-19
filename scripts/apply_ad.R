suppressMessages({
    library(dplyr)
    library(kSamples)
    library(purrr)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_ad <- function(sce, pars, ds_only = TRUE) {
    t <- system.time({
        if (!ds_only) 
            assays(sce, pars$assay) <- switch(pars$assay,
                logcounts = logNormCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
        res <- tryCatch(
            error = function(e) e, 
            run_ad(sce, pars$assay, pars$var))
    })[[3]]
    list(rt = t, tbl = res)
}

run_ad <- function(sce, assay, var) {
    form <- as.formula(paste("y ~", var))
    kids <- levels(sce$cluster_id)
    cs_by_k <- split(colnames(sce), sce$cluster_id)
    res <- lapply(kids, function(k) {
        sub <- sce[, cs_by_k[[k]]]
        sub <- sub[rowSums(assay(sub) > 0) >= 10, ]
        apply(assay(sub, assay), 1, function(y)
            ad.test(form, data = data.frame(y, colData(sub)))) %>% 
            map(function(u) u$ad[5]) %>% 
            unlist %>% data.frame(
                gene = names(.), 
                cluster_id = k,
                p_val = .,
                p_adj.loc = p.adjust(.),
                row.names = NULL,
                stringsAsFactors = FALSE)
    })
    df <- bind_rows(res)
    df$p_adj.glb <- p.adjust(df$p_val)
    return(df)
}

