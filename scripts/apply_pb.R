apply_pb <- function(sce, pars) {
    
    # specify design & contrast matrix
    ei <- metadata(sce)$experiment_info
    design <- model.matrix(~ 0 + ei$group_id)
    dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
    contrast <- makeContrasts("B-A", levels = design)
    
    # get aggregateData() parameters
    args <- as.list(args("aggregateData"))
    pb_pars <- pars[names(pars) %in% names(args)]
    pb_pars <- c(list(sce), pb_pars)
    
    # run & time method
    t1 <- system.time(pb <- do.call(aggregateData, pb_pars))[[3]]
    
    # get runDS() parameters
    args <- as.list(args("pbDS"))
    ds_pars <- pars[names(pars) %in% names(args)]
    ds_pars <- c(list(sce, pb, design, contrast, verbose = FALSE), ds_pars)
    
    # run & time method
    t2 <- system.time(ds <- do.call(pbDS, ds_pars))[[3]]
    res <- bind_rows(ds$table[[1]])
    
    # return results
    list(rt = c(t1, t2), res = res)
}

