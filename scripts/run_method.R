# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(limma)
        library(muscat)
        library(purrr)
        library(SingleCellExperiment)
    }))

# load data
sce <- readRDS(snakemake@input$sim_data)

# load parameters
pars <- snakemake@input$method_pars
pars <- read.csv(pars, row.names = 1, stringsAsFactors = FALSE)
pars <- as.list(pars[snakemake@wildcards$method_id, ])

# load configuration
L <- jsonlite::fromJSON(snakemake@input$config)

# source method function
source(fun <- snakemake@input$apply_fun)
fun <- gsub("(.R)", "", basename(fun))

# split cells by cluster-sample
cells_by_ks <- muscat:::.split_cells(sce)

g <- snakemake@wildcards$g
c <- snakemake@wildcards$c
s <- snakemake@wildcards$s
i <- snakemake@wildcards$i_run
set.seed(L$seed + as.numeric(i))

# get genes to use
genes_use <- switch(g, 
    all = seq_len(nrow(sce)), 
    sample(nrow(sce), g))
# subset samples to use
sids <- levels(sce$sample_id)
m <- match(sids, sce$sample_id)
gids <- sce$group_id[m]
samples_use <- switch(s,
    all = sids,
    sapply(split(sids, gids), sample, s))
cells_by_ks <- lapply(cells_by_ks, "[", samples_use)
# get cells to use
cells_use <- switch(c, 
    all = unlist(cells_by_ks), 
    map_depth(cells_by_ks, 2, sample, c) %>% unlist)

# subset, apply method & write results to .rds
sub <- sce[genes_use, cells_use]
res <- get(fun)(sub, pars)
saveRDS(res, snakemake@output$res)


