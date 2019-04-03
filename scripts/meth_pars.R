# load packages
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

# pbDS() parameters
pb_pars <- expand.grid(
    stringsAsFactors = FALSE,
    assay = c("counts", "normcounts", "logcounts", "cpm"),
    fun = c("sum", "mean", "median"),
    scale = c(TRUE, FALSE),
    method = c("edgeR", "limma-voom", "limma-trend")) %>% 
    filter(
        !(assay != "cpm" & scale), # scale only w/ CPM
        !(assay == "counts" & fun == "mean"), # mean not w/ raw counts
        !(!(assay == "cpm" & scale) & fun == "median"), # median only w/ scaledCPM
        !(assay != "counts" & method %in% c("edgeR", "limma-voom"))) # edgeR & limma-voom only w/ counts

pb_pars$id <- with(pb_pars, sprintf(sprintf("%s(%s%s).%s", 
    fun, c("", "scale")[as.numeric(scale)+1], assay, method)))

# mixed-model parameters -------------------------------------------------------
mm_pars <- expand.grid(
    stringsAsFactors = FALSE,
    method = c("dream", "vst"),
    covs = c("", "dr"))

mm_pars$id <- with(mm_pars, {
    sep <- c("", ".")[as.numeric(covs != "")+1]
    sprintf("mm-%s%s%s", method, sep, covs)
})

# MAST parameters --------------------------------------------------------------
mast_pars <- expand.grid(
    covs = c("", "dr"))

mast_pars$id <- with(mast_pars, {
    sep <- c("", ".")[as.numeric(covs != "")+1]
    sprintf("MAST%s%s", sep, covs)
})

# scDD & AD parameters ---------------------------------------------------------
scdd_pars <- data.frame(id = "scDD")
ad_pars <- data.frame(id = "AD")

# write parameters to .csv
for (i in c("pb", "mm", "ad", "mast", "scdd")) {
    j <- paste(i, "pars", sep = "_")
    fn <- file.path("metadata", paste0(j, ".csv"))
    write.csv(get(j), fn, row.names = get(j)$id)
}

# write method IDs to .rds
method_ids <- lapply(c("pb", "mm", "mast"), function(i) {
    ids <- get(paste(i, "pars", sep = "_"))$id
    setNames(rep(i, length(ids)), ids) 
}) %>% unlist() %>% 
    data.frame(
        type = ., id = names(.), 
        stringsAsFactors = FALSE) %>% 
    bind_rows(
        data.frame(
            row.names = "scDD", 
            type = "scdd", id = "scDD",
            stringsAsFactors = FALSE),
        data.frame(
            row.names = "AD",
            type = "ad", id = "AD",
            stringsAsFactors = FALSE)) %>% 
    write.csv(file.path("metadata", "method_ids.csv"))

