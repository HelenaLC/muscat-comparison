suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(kSamples)
        library(purrr)
        library(scater)
        library(sctransform)
        library(SingleCellExperiment)
    })
)

apply_ad <- function(sce, pars, ds_only = TRUE) {
    # run & time method
    t <- system.time({
        if (!ds_only) {
            suppressWarnings(suppressMessages(
                assays(sce)[[pars$assay]] <- switch(pars$assay,
                    logcounts = logcounts(normalize(sce)),
                    vstcounts = vst(counts(sce), show_progress = FALSE)$y)))
        }
        res <- tryCatch(
            error = function(e) e, 
            run_ad(sce, pars$assay, pars$var))
    })[[3]]
    
    # return results
    list(rt = t, tbl = res)
}

run_ad <- function(sce, assay, var) {
    form <- as.formula(paste("y ~", var))
    kids <- levels(sce$cluster_id)
    cs_by_k <- split(colnames(sce), sce$cluster_id)
    lapply(kids, function(k) {
        sub <- sce[, cs_by_k[[k]]]
        sub <- sub[rowSums(assay(sub) > 0) > 10, ]
        apply(assays(sub)[[assay]], 1, function(y)
            ad.test(form, data = data.frame(y, colData(sub)))) %>% 
            map(function(u) u$ad[5]) %>% 
            unlist %>% data.frame(stringsAsFactors = FALSE,
                gene = names(.), cluster_id = k, p_val = .) %>% 
            dplyr::mutate(p_adj.loc = p.adjust(p_val))
    }) %>% bind_rows %>% dplyr::mutate(p_adj.glb = p.adjust(p_val))
}

