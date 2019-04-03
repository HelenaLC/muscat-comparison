apply_mast <- function(sce, pars) {
    
    # get method parameters
    args <- as.list(args("run_mast"))
    pars <- pars[names(pars) %in% names(args)]
    
    # compute logCPM
    assays(sce)$logcpm <- log2(assays(sce)$cpm + 1)
    
    # run & time method
    t <- system.time({
        if (pars$covs == "dr") {
            sce$dr <- colMeans(assays(sce)$counts > 0)
        } else {
            pars$covs <- NULL
        }
        res <- tryCatch(error = function(e) NULL, 
            do.call(run_mast, c(list(sce), pars)))
    })[[3]]
    
    # return results
    list(rt = t, res = res)
}

# ------------------------------------------------------------------------------
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(edgeR)
        library(MAST)
        library(SummarizedExperiment)
    })
)

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
        p_adj.loc <- stats::p.adjust(p_val)
        data.frame(
            gene = rownames(y),
            cluster_id = k,
            p_val, p_adj.loc,
            row.names = NULL,
            stringsAsFactors = FALSE)
    }) %>% bind_rows %>% mutate(p_adj.glb = p.adjust(p_val))
}