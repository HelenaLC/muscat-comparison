suppressMessages({
    library(dplyr)
    library(MAST)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_mast <- function(sce, pars, ds_only = TRUE) {
    # run & time method
    t <- system.time({
        if (!ds_only) {
            suppressWarnings(
                assays(sce)$lCount <- switch(pars$assay, 
                    logcounts = logcounts(normalize(sce)),
                    vstresiduals = vst(counts(sim), show_progress = FALSE)$y))
        } else {
            assays(sce)$lCount <- assays(sce)[[pars$assay]]
        }
        res <- tryCatch(
            error = function(e) e, 
            run_mast(sce))
    })[[3]]
    
    # return results
    list(rt = t, tbl = res)
}

run_mast <- function(sce) {
    colData(sce)$wellKey <- colnames(sce)
    rowData(sce)$primerid <- rownames(sce)
    
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    kids <- levels(sce$cluster_id)
    
    lapply(kids, function(k) {
        y <- sce[, cells_by_k[[k]]]
        sca <- SceToSingleCellAssay(y, check_sanity = FALSE)
        fit <- suppressMessages(zlm(~ group_id, sca))
        lrt <- suppressMessages(lrTest(fit, "group_id"))
        p_val <- lrt[, "hurdle", "Pr(>Chisq)"]
        p_adj.loc <- p.adjust(p_val)
        data.frame(
            gene = rownames(y),
            cluster_id = k,
            p_val, p_adj.loc,
            row.names = NULL,
            stringsAsFactors = FALSE)
    }) %>% bind_rows %>% dplyr::mutate(p_adj.glb = p.adjust(p_val))
}
