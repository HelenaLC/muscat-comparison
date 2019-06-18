config <- yaml::read_yaml("config.yaml")
ids <- purrr::set_names(c("pb", "ad", "scdd", "mast", "mm"))

# aggregation-based ------------------------------------------------------------
pb <- dplyr::bind_rows(
    expand.grid(
        stringsAsFactors = FALSE,
        assay = "counts", fun = "sum", scale = FALSE, 
        method = c("edgeR", "limma-voom", "limma-trend")
    ),
    expand.grid(
        stringsAsFactors = FALSE, scale = FALSE,
        assay = c("normcounts", "cpm", "logcpm"),
        fun = "sum", method = "limma-trend"
    ),
    expand.grid(
        stringsAsFactors = FALSE, scale = FALSE,
        assay = c("normcounts", "logcounts", "vstcounts"),
        fun = "mean", method = "limma-trend"
    ),
    data.frame(
        stringsAsFactors = FALSE,
        assay = "cpm", fun = "median", scale = TRUE, method = "limma-trend"
    )    
)
pb$id <- with(pb, sprintf("%s.%s(%s%s)", 
    method, fun, ifelse(scale, "scale", ""), assay))

# mixed-models -----------------------------------------------------------------
mm <- dplyr::bind_rows(
    expand.grid(
        stringsAsFactors = FALSE,
        method = "dream",
        vst = "",
        covs = c("dr", ""),
        ddf = "Kenward-Roger"),
    expand.grid(
        stringsAsFactors = FALSE,
        method = "vst",
        vst = c("sctransform", "DESeq2"),
        covs = c("dr", ""),
        ddf = "Kenward-Roger"))
mm$id <- with(mm, {
    sep1 <- ifelse(vst != "", ".", "")
    sep2 <- ifelse(covs != "", "_", "")
    sprintf("MM-%s%s%s%s%s", method, sep1, vst, sep2, covs)
})

# Anderson-Darling -------------------------------------------------------------
ad <-  expand.grid(
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstcounts"),
    var = c("sample_id", "group_id"))
ad$id <- with(ad, sprintf("AD-%s.%s",
    gsub("(.).*", "\\1id", var), assay))

# scDD -------------------------------------------------------------------------
scdd <- data.frame(
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstcounts"))
scdd$id <- with(scdd, sprintf("scDD.%s", assay))

# MAST -------------------------------------------------------------------------
mast <- expand.grid(
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE,
    assay = c("logcounts", "vstcounts"),
    covs = c("dr", ""))
mast$id <- with(mast, sprintf("MAST.%s%s%s",
    assay, ifelse(covs == "", "", "_"), covs))

# write method IDs to .csv -----------------------------------------------------
for (id in ids) {
    ps <- get(id)
    assign(id, split(ps, ps$id))
}

pars <- sapply(ids, get)
typs <- rep.int(ids, sapply(pars, length))
pars <- purrr::flatten(pars)
mids <- names(pars)

mids_df <- data.frame(
    row.names = NULL,
    stringsAsFactors = FALSE,
    id = mids, type = typs)

write.csv(mids_df, config$mids)

# write method parameters to .json ---------------------------------------------
fns <- paste0(mids, ".json")
fns <- file.path(config$meth_pars, fns)
names(fns) <- mids

for (id in mids) {
    new <- as.list(pars[[id]])
    if (file.exists(fns[id])) {
        old <- jsonlite::fromJSON(fns[id])
        if (!identical(old, new))
            write(jsonlite::toJSON(new), fns[id])
    } else {
        write(jsonlite::toJSON(new), fns[id])
    }
}