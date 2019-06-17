suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(muscat)
        library(scater)
        library(sctransform)
        library(SingleCellExperiment)
    })
)
source(snakemake@config$utils)

sce <- readRDS(snakemake@input$sce)
sim_pars <- jsonlite::fromJSON(snakemake@input$sim_pars)

i <- as.numeric(snakemake@wildcards$i)
set.seed(sim_pars$seed + i)

kids <- levels(sce$cluster_id)
sids <- levels(sce$sample_id)

sce <- .filter_sce(sce, 
    kids = sample(kids, sim_pars$nk),
    sids = sample(sids, sim_pars$ns))

sim <- simData(sce, nrow(sce), sim_pars$nc, p_dd = sim_pars$p_dd)
sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
sim <- sim[sample(nrow(sim), sim_pars$ng), ]

gi <- metadata(sim)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim))
metadata(sim)$gene_info <- gi

sim$dr <- colMeans(counts(sim) > 0)
assays(sim)$cpm <- calculateCPM(sim)
assays(sim)$logcpm <- log2(cpm(sim)+1)
assays(sim)$vstcounts <- suppressWarnings(vst(counts(sim), show_progress = FALSE)$y)
sim <- suppressWarnings(normalize(normalize(sim, return_log = FALSE)))

saveRDS(sim, snakemake@output$sim)