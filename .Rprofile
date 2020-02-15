options(conflicts.policy = list(warn = FALSE))
data.table::setDTthreads(threads = 1)

# source utilities
config <- yaml::read_yaml("config.yaml")
source(config$utils)

# wrapper to deprase snakemake wildcards
.get_wcs <- function(wcs) {
    ss <- strsplit(wcs, ",")[[1]]
    ss <- sapply(ss, strsplit, "=")
    keys <- sapply(ss, .subset, 1)
    vals <- sapply(ss, .subset, 2)
	wcs <- as.list(vals)
	names(wcs) <- keys
    return(wcs)
}

# get commandline arguments
args <- R.utils::commandArgs(
	trailingOnly = TRUE, 
	asValues = TRUE)

if (!is.null(args$wcs))
	wcs <- .get_wcs(args$wcs)

for (i in seq_along(args))
    if (any(grepl(";", args[[i]])))
        args[[i]] <- unlist(strsplit(args[[i]], ";"))

# print input arguments for logs
if (length(args) > 0) print(args)
if (exists("wcs")) print(wcs)
