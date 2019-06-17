suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    # specify design & contrast matrix
    ei <- metadata(sce)$experiment_info
    design <- model.matrix(~ 0 + ei$group_id)
    dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
    contrast <- makeContrasts("B-A", levels = design)
    
    # get aggregateData() parameters
    args <- as.list(args("aggregateData"))
    agg_pars <- pars[names(pars) %in% names(args)]
    
    # run & time method
    t1 <- system.time({
        if (!ds_only) {
            a <- agg_pars$assay
            suppressWarnings(suppressMessages(
                assays(sce)[[a]] <- switch(a, 
                    counts = counts(sce),
                    normcounts = normcounts(normalize(sce, return_log = FALSE)),
                    logcounts = logcounts(normalize(sce)),
                    vstcounts = vst(counts(sce), show_progress = FALSE)$y,
                    cpm = calculateCPM(counts(sce)),
                    logcpm = log2(calculateCPM(counts(sce)) + 1))))
        }
        pb <- do.call(aggregateData, c(list(sce), agg_pars))
    })[[3]]
    
    # get pbDS() parameters
    args <- as.list(args("pbDS"))
    ds_pars <- pars[names(pars) %in% names(args)]
    ds_pars <- c(list(sce, pb, design, contrast, verbose = FALSE), ds_pars)
    
    # run & time method
    t2 <- system.time(
        suppressWarnings(
            res <- tryCatch(
                error = function(e) NULL, 
                do.call(pbDS, ds_pars))))[[3]]
    
    # return results
    list(rt = c(t1, t2), tbl = dplyr::bind_rows(res$table[[1]]))
}