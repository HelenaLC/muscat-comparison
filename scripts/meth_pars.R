config <- yaml::read_yaml("config.yaml")
method_ids <- c("pb", "ad", "mast", "scdd")

# load packages
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

# pbDS() parameters
pb_pars <- expand.grid(
    stringsAsFactors = FALSE,
    assay = c("counts", "normcounts", "logcounts", "vstcounts", "cpm"),
    fun = c("sum", "mean", "median"),
    scale = c(TRUE, FALSE),
    method = c("edgeR", "limma-voom", "limma-trend")) %>% 
    filter(
        !(assay != "cpm" & scale), # scale only w/ CPM
        !(assay == "counts" & fun == "mean"), # mean not w/ raw counts
        !(!(assay == "cpm" & scale) & fun == "median"), # median only w/ scaledCPM
        !(assay != "counts" & method %in% c("edgeR", "limma-voom")), # edgeR & limma-voom only w/ counts
        !(fun == "sum" & method == "limma-trend" & assay %in% c("logcounts", "vstcounts")))

pb_pars$id <- with(pb_pars, sprintf(sprintf("%s.%s(%s%s)", 
    method, fun, c("", "scale")[as.numeric(scale)+1], assay)))

# mixed-model parameters -------------------------------------------------------
mm_pars <- rbind(
    expand.grid(
        stringsAsFactors = FALSE,
        method = "vst",
        vst = c("sctransform", "DESeq2"),
        covs = c("dr", ""),
        ddf = "Kenward-Roger"),
    expand.grid(
        stringsAsFactors = FALSE,
        method = "dream",
        vst = "",
        covs = c("dr", ""),
        ddf = "Kenward-Roger"))

mm_pars$id <- with(mm_pars, {
    sep1 <- c("", ".")[as.numeric(vst != "")+1]
    sep2 <- c("", "_")[as.numeric(covs != "")+1]
    sprintf("MM-%s%s%s%s%s", method, sep1, vst, sep2, covs)
})

# MAST parameters --------------------------------------------------------------
mast_pars <- expand.grid(
    stringsAsFactors = FALSE,
    assay = c("logcpm", "logcounts"),
    covs = c("", "dr"))

mast_pars$id <- with(mast_pars, {
    sep <- c("", "_")[as.numeric(covs != "")+1]
    sprintf("MAST.%s%s%s", assay, sep, covs)
})

# scDD parameters --------------------------------------------------------------
scdd_pars <- data.frame(
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstcounts"))
scdd_pars$id <- with(scdd_pars, paste("scDD", assay, sep = "."))

# AD parameters ----------------------------------------------------------------
ad_pars <- data.frame(
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstcounts"))
ad_pars$id <- with(ad_pars, paste("AD", assay, sep = "."))

# write parameters to .csv
for (i in paste(method_ids, "pars", sep = "_")) {
    pars_new <- get(i)
    rownames(pars_new) <- pars_new$id
    fn <- paste(i, "csv", sep = ".")
    fn <- file.path(config$meth_pars, fn)
    if (file.exists(fn)) {
        pars_old <- read.csv(fn, row.names = 1, stringsAsFactors = FALSE)
        if (!isTRUE(all_equal(pars_new, pars_old)))   # only write new file if
            write.csv(pars_new, fn, row.names = TRUE) # something has changed
    } else {
        write.csv(pars_new, fn, row.names = TRUE)
    }
}

# write method annotations to .csv
lapply(method_ids, function(id) {
    pars <- get(paste(id, "pars", sep = "_"))
    typ <- switch(id, mm = "cell", pb = "sample", "group")
    inp <- switch(id, 
        mm = vapply(pars$method, function(method)
            switch(method,
                dream = "logcounts", 
                vst = "vstcounts"),
            character(1)),
        pars$assay)
    agg <- switch(id, pb = pars$fun, "none")
    log <- switch(id, pb = grepl("log", inp), TRUE)
    data.frame(
        id = as.character(pars$id), 
        typ, inp, log, agg, 
        stringsAsFactors = FALSE)
}) %>% bind_rows %>% 
    mutate_at("typ", paste, "level", sep = "-") %>% 
    mutate_at("typ", factor, levels = paste(c(
        "cell", "sample", "group"), "level", sep = "-")) %>% 
    mutate_at("agg", factor, levels = c("sum", "mean", "median", "none")) %>% 
    mutate_at("log", factor, levels = c("TRUE", "FALSE")) %>% 
    saveRDS(config$anno_tbl)

# write method IDs to .rds
lapply(method_ids, function(i) {
    ids <- get(paste(i, "pars", sep = "_"))$id
    setNames(rep(i, length(ids)), ids) 
}) %>% unlist() %>% data.frame(
    type = ., id = names(.), 
    stringsAsFactors = FALSE) %>% 
    write.csv(config$method_ids)

