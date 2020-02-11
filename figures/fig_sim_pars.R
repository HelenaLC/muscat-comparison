suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(muscat)
    library(purrr)
    library(scater)
    library(scran)
    library(SingleCellExperiment)
})

set.seed(160185)

# sce <- readRDS(file.path(config$raw_data, "ref_kang.rds"))
# sce <- .filter_sce(sce, 
#     kids = c("B cells", "CD14+ Monocytes", "CD4 T cells"), 
#     sids = c("ctrl1015", "ctrl1256"))

sim_pars <- dplyr::bind_rows(
    expand.grid(
        p_type = c(0, 0.02, 0.05),
        p_dd = list(diag(6)[1, ]),
        probs = list(list(NULL, NULL, c(1, 0)))
    ),
    expand.grid(
        p_type = 0,
        p_dd = list(
            c(0.98, 0, 0.02, rep(0, 3)),
            c(0.95, 0, 0.05, rep(0, 3)), 
            c(0.90, 0, 0.10, rep(0, 3)))
    ),
    expand.grid(
        p_type = 0,
        p_dd = list(c(0.98, 0, 0.02, rep(0, 3))),
        rel_lfc = list((c(1, 1, 1)), c(0.5, 1, 1.5), c(0, 1, 2))
    )
)

sim <- lapply(seq_len(nrow(sim_pars)), function(i) {
    u <- as.list(sim_pars[i, ])
    u <- purrr::map(u, 1)
    u <- do.call(simData, c(u, list(x = sce, 
        ns = 3, nk = 3, nc = 2*2*3*100)))
    u <- u[sample(nrow(u), 4e3), ]
    u <- logNormCounts(u)
})

sim <- lapply(sim, function(u) {
    u <- runPCA(u, ncomponents = 20)
    u <- runTSNE(u, use_dimred = "PCA", n_dimred = 20)
})

labs <- c(
    sprintf("%s%% DE, 0%% DS", with(sim_pars[1:3, ], p_type * 100)),
    sprintf("%s%% DS, 0%% DE", with(sim_pars[4:6, ], unlist(map(p_dd, 3)) * 100)),
    sprintf("logFC %s, 5%% DS", sapply(map(sim_pars[7:9, "rel_lfc"], `*`, 2), paste, collapse = "/")))

df <- map(sim, function(u) {
    dr <- reducedDim(u, "TSNE")
    df <- data.frame(dr, colData(u)) %>% 
        mutate_if(is.factor, as.character)
}) %>% bind_rows(.id = "i") %>% 
    dplyr::mutate(id = factor(paste(cluster_id, group_id))) %>% 
    mutate_at("i", factor, labels = labs)

cols <- CATALYST:::.cluster_cols[
    seq_len(nlevels(df$id))] %>% 
    set_names(levels(df$id))

p <- ggplot(df, aes(x = X1, y = X2, col = id)) +
    facet_wrap(~ i, ncol = 3, scales = "free") + 
    geom_point_rast(size = 0.2, alpha = 0.4, raster.dpi = 50) + 
    scale_color_manual(values = cols) +
    .prettify("void") + theme(
        legend.position = "none",
        panel.spacing = unit(2, "mm"),
        panel.border = element_rect(),
        axis.text = element_blank(),
        strip.text = element_text(face = "bold", size = 4, hjust = 0))

ggsave(file.path("figures", "sim_pars.pdf"), p,
    width = 6, height = 6.5, units = "cm",
    dpi = 300, useDingbats = FALSE)










