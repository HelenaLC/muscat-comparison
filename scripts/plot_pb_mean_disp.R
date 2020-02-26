suppressMessages({
    library(dplyr)
    library(edgeR)
    library(ggplot2)
    library(ggrastr)
    library(muscat)
    library(SingleCellExperiment)
})

set.seed(1994)
sce <- readRDS(args$sce)

nk <- ns <- 4
sim <- simData(sce, nk = nk, ns = ns,
    nc = 2*nk*ns*200, ng = nrow(sce))

pbs <- list(
    reference = aggregateData(sce),
    simulation = aggregateData(sim))

assays(pbs$reference) <- assays(pbs$reference)[metadata(sim)$ref_kids]
assayNames(pbs$reference) <- assayNames(pbs$simulation)

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

df <- dplyr::filter(res, 
    disp_tag > quantile(disp_tag, 0.01),
    disp_tag < quantile(disp_tag, 0.99))

p <- ggplot(df, aes(x = mean, col = id)) +
    facet_wrap(~ cluster_id, nrow = 1) + 
    scale_x_continuous(
        limits = c(0, 18), 
        breaks = seq(0, 18, 6), 
        expand = c(0.05, 0)) +
    scale_y_continuous(
        trans = "log10", 
        limits = c(0.01, 0.3), 
        expand = c(0.05, 0)) +
    geom_point_rast(aes(y = disp_tag), size = 0.2, alpha = 0.02) + 
    geom_point_rast(aes(y = disp_fit), size = 0.1, alpha = 0.1) +
    scale_color_manual(NULL, values = c("royalblue", "tomato")) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "mean logCPM", y = "dispersion") +
    .prettify("bw", legend.position = "bottom")

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 6, units = "cm", 
    dpi = 300, useDingbats = FALSE)

