suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(MAST)
        library(scater)
        library(SingleCellExperiment)
    })
)

apply_mast <- function(sce, pars, ds_only = TRUE) {
    # run & time method
    t <- system.time({
        if (pars$covs == "dr") {
            sce$dr <- colMeans(assays(sce)$counts > 0)
        } else {
            pars$covs <- NULL
        }
        if (!ds_only) {
            assays(sce)$lCount <- switch(pars$assay, 
                logcpm = log2(calculateCPM(sce) + 1),
                logcounts = suppressWarnings(logcounts(normalize(sce))))
        } else {
            assays(sce)$lCount <- assays(sce)[[pars$assay]]
        }
        res <- tryCatch(
            error = function(e) NULL, 
            run_mast(sce, pars$covs))
    })[[3]]
    
    # return results
    list(rt = t, res = res)
}

run_mast <- function(sce, covs = NULL) {
    formula <- ~ group_id
    if (!is.null(covs)) {
        formula <- paste(as.character(formula), collapse = "")
        formula <- paste(formula, covs, sep = "+")
        formula <- as.formula(formula)
    }
    
    colData(sce)$wellKey <- colnames(sce)
    rowData(sce)$primerid <- rownames(sce)
    
    cells_by_k <- split(colnames(sce), sce$cluster_id)
    kids <- levels(sce$cluster_id)
    
    lapply(kids, function(k) {
        y <- sce[, cells_by_k[[k]]]
        sca <- SceToSingleCellAssay(y, check_sanity = FALSE)
        fit <- suppressMessages(zlm(formula, sca))
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
