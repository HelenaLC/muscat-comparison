library(dplyr)
library(edgeR)
library(ggplot2)

source("scripts/utils.R")
sce <- readRDS(snakemake@input$sce)
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, dim = c(4, 2))
sce <- .filter_sce(sce, dimnames(n_cells))
sim <- simData(sce, nrow(sce), 2400, probs = list(NULL, NULL, c(1, 0)))

pb1 <- aggregateData(sce)
pb2 <- aggregateData(sim)
pbs <- setNames(list(pb1, pb2), c("Kang", "sim"))

res <- lapply(names(pbs), function(id) {
    u <- pbs[[id]]
    lapply(assayNames(u), function(k) {
        y <- assays(u)[[k]]
        y <- DGEList(y)
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

p <- ggplot(res, aes(x = mean, col = id)) +
    facet_wrap(~ cluster_id, nrow = 1) + 
    geom_point(aes(y = disp_tag), size = 0.01, alpha = 0.01) + 
    geom_point(aes(y = disp_fit), size = 0.02, alpha = 0.1) +
    scale_color_manual(values = c("navy", "maroon"), guide = FALSE) +
    labs(x = "mean logCPM", y = "dispersion") +
    scale_x_continuous(limits = c(1, 15), breaks = seq(2, 14, 4), expand = c(0, 0.2)) +
    scale_y_log10(limits = c(1e-4, 1), breaks = 10^seq(-4, 0, 1), expand = c(0, 0.3)) + 
    prettify() + theme(aspect.ratio = 3/2)
ggsave("figures/sala/pb_mean_disp.pdf", p, width = 15, height = 6.5, unit = "cm", dpi = 300)

   
    