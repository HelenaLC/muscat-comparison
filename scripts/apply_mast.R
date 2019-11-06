suppressMessages({
    library(dplyr)
    library(MAST)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_mast <- function(sce, pars, ds_only = TRUE) {
    t <- system.time({
        if (ds_only) {
            assay(sce, "lCount") <- assay(sce, pars$assay)
        } else {
            assay(sce, "lCount") <- switch(pars$assay, 
                logcounts = logNormCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sim), show_progress = FALSE)$y)
        }
        res <- tryCatch(run_mast(sce), error = function(e) e)
    })[[3]]
    list(rt = t, tbl = res)
}

run_mast <- function(sce) {
    colData(sce)$wellKey <- colnames(sce)
    rowData(sce)$primerid <- rownames(sce)
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    kids <- levels(sce$cluster_id)
    res <- lapply(kids, function(k) {
        y <- sce[, cells_by_k[[k]]]
        sca <- SceToSingleCellAssay(y, check_sanity = FALSE)
        fit <- suppressMessages(zlm(~ group_id, sca))
        lrt <- suppressMessages(lrTest(fit, "group_id"))
        p_val <- lrt[, "hurdle", "Pr(>Chisq)"]
        data.frame(
            gene = rownames(y), 
            cluster_id = k, 
            p_val, 
            p_adj.loc = p.adjust(p_val), 
            row.names = NULL, 
            stringsAsFactors = FALSE)
    }) 
    df <- bind_rows(res)
    df$p_adj.glb <- p.adjust(df$p_val)
    return(df)
}