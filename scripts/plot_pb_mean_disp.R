source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(dplyr)
    library(edgeR)
    library(ggplot2)
    library(muscat)
    library(SingleCellExperiment)
})

set.seed(1994)

sce <- readRDS(snakemake@input$sce)

kids <- levels(sce$cluster_id)
sids <- levels(sce$sample_id)

sce <- .filter_sce(sce, 
    kids = sample(kids, (nk <- 4)), 
    sids = sample(sids, (ns <- 3)))

sim <- simData(sce, nrow(sce), 2*nk*ns*200)

pbs <- list(
    reference = aggregateData(sce),
    simulation = aggregateData(sim))

res <- lapply(names(pbs), function(id) {
    u <- pbs[[id]]
    lapply(assayNames(u), function(k) {
        y <- assays(u)[[k]]
        y <- DGEList(y, remove.zeros = TRUE)
        y <- calcNormFactors(y)
        y <- estimateDisp(y)
        data.frame(id,
            cluster_id = k,
            mean = y$AveLogCPM, 
            disp_tag = y$tagwise.dispersion,
            disp_fit = y$trended.dispersion,
            stringsAsFactors = FALSE)
    }) %>% bind_rows
}) %>% bind_rows

kids <- paste0("cluster", seq_len(nk))
kids <- kids[as.numeric(factor(res$cluster_id))]
res$cluster_id <- kids

res <- dplyr::filter(res, 
    disp_tag > quantile(disp_tag, 0.01),
    disp_tag < quantile(disp_tag, 0.99))

sub <- res[sample(nrow(res), 1e4), ]

p <- ggplot(sub, aes(x = mean, col = id)) +
    facet_wrap(~ cluster_id, nrow = 1) + 
    geom_point(aes(y = disp_tag), size = 0.2, alpha = 0.02) + 
    geom_point(aes(y = disp_fit), size = 0.1, alpha = 0.1) +
    scale_color_manual(NULL, values = c("royalblue", "tomato")) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "mean logCPM", y = "dispersion") + scale_y_log10() + 
    .prettify("bw", legend.position = "bottom")

saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    width = 15, height = 6, units = "cm", 
    dpi = 300, useDingbats = FALSE)