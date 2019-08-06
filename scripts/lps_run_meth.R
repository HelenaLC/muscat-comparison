options(conflicts.policy = list(warn = FALSE))

# load data & method parameters
sce <- readRDS(file.path("MAGL", "output", "MAGL-SCE.rds"))
meth_pars <- as.list(jsonlite::fromJSON(snakemake@input$meth_pars))

# get method wrapper 'apply_xx()'
source(fun <- snakemake@input$fun)
fun <- gsub("(.R)", "", basename(fun))

# run method & write results to .rds
res <- get(fun)(sce, meth_pars)
saveRDS(res, snakemake@output$res)