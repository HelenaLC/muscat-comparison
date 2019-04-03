# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(muscat)
        library(scater)
        library(SingleCellExperiment)
    }))

# source helpers
source(snakemake@config$utils)

# load data
sce <- readRDS(snakemake@input$data)

# get simulation parameters
sim_pars <- jsonlite::fromJSON(snakemake@input$sim_pars)

# default to 3 clusters & 3 samples
if (!"nk" %in% names(sim_pars)) sim_pars$nk <- 3
if (!"ns" %in% names(sim_pars)) sim_pars$ns <- 3

# filter clusters & samples
n_cells <- as.matrix(table(sce$cluster_id, sce$sample_id))
n_cells <- .filter_matrix(n_cells, dim = c(sim_pars$nk, sim_pars$ns))

# set seed (input + rep. nb.)
i <- snakemake@wildcards$i_sim
set.seed(sim_pars$seed + as.numeric(i))

# subset SCE
sce <- .filter_sce(sce, dimnames(n_cells))
sim_pars <- sim_pars[names(sim_pars) %in% 
        names(as.list(args("simData")))]

# simulate, normalize, and compute CPM
sim <- do.call(simData, c(
    list(x = sce, n_genes = nrow(sce)), 
    sim_pars[names(sim_pars) != "n_genes"]))
suppressWarnings(sim <- normalize(normalize(sim, return_log = FALSE)))
assays(sim)$cpm <- calculateCPM(sim)
assays(sim)$logcpm <- log2(assays(sim)$cpm + 1)

# downsample to n_genes & write to .rds
gs <- sample(rownames(sim), sim_pars$n_genes)
gi <- metadata(sim)$gene_info
metadata(sim)$gene_info <- gi %>% dplyr::filter(gene %in% gs)
sim <- sim[gs, ]
saveRDS(sim, snakemake@output$sim_data)
