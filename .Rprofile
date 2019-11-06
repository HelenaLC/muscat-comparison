readRenviron(".Renviron")
.libPaths(c(Sys.getenv("R_LIBS")))
options(conflicts.policy = list(warn = FALSE))
if (exists("snakemake"))
	source(snakemake@config$utils)

