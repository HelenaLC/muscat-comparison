suppressMessages({
    library(BiocParallel)
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
            assay(sce, pars$assay) <- switch(pars$assay,
                logcounts = normalizeCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
        if (!is.matrix(a <- assay(sce, pars$assay))) 
            assay(sce, pars$assay) <- as.matrix(a)
        if (is.null(n <- pars$n_threads)) n <- 1
        res <- tryCatch(
            error = function(e) e, 
            run_ad(sce, pars$assay, pars$var, n))
    })[[3]]
    list(rt = t, tbl = res)
}

run_ad <- function(sce, assay, var, n_threads) {
    form <- as.formula(paste("y ~", var))
    kids <- levels(sce$cluster_id)
    cs_by_k <- split(colnames(sce), sce$cluster_id)
    lapply(kids, function(k) {
        y <- sce[, cs_by_k[[k]]]
        y <- sub[rowSums(assay(sub) > 0) >= 10, ]
        bplapply(seq_len(nrow(y)), function(g) {
            ad.test(form, data = data.frame(
                y = assay(y, assay)[g, ], 
                colData(y)))$ad[5]
        }, BPPARAM = MulticoreParam(n_threads)) %>% 
            unlist %>% data.frame(
                gene = rownames(y), cluster_id = k,
                p_val = ., p_adj.loc = p.adjust(.),
                row.names = NULL, stringsAsFactors = FALSE)
    }) %>% bind_rows %>% mutate(p_adj.glb = p.adjust(p_val))
}