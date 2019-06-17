options(conflicts.policy = list(warn = FALSE))

source(snakemake@config$utils)

suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(muscat)
    library(purrr)
    library(scater)
    library(scran)
    library(SingleCellExperiment)
})

set.seed(1601)

sce <- readRDS("data/raw_data/kang.rds")

sce <- .filter_sce(sce, 
    kids = c("B cells", "CD14+ Monocytes", "CD4 T cells"), 
    sids = c("ctrl1015", "ctrl1488"))

sim_pars <- dplyr::bind_rows(
    expand.grid(
        p_type = c(0, 0.01, 0.05),
        p_dd = list(diag(6)[1, ]),
        probs = list(list(NULL, NULL, c(1, 0)))
    ),
    expand.grid(
        p_type = 0,
        p_dd = list(diag(6)[1, ], 
            c(0.99, 0, 0.01, rep(0, 3)), 
            c(0.95, 0, 0.05, rep(0, 3)))
    ),
    expand.grid(
        p_type = 0,
        p_dd = list(c(0.95, 0, 0.05, rep(0, 3))),
        rel_lfc = list((c(1, 1, 1)), c(0.8, 1, 1.2), c(0.5, 1, 1.5))
    )
)

sim <- lapply(row(sim_pars)[, 1], function(i) {
    u <- as.list(sim_pars[i, ])
    u <- purrr::map(u, 1)
    do.call(simData, c(u, list(x = sce, 
        n_genes = nrow(sce), n_cells = 2*2*3*100)))
})

sim2 <- lapply(sim, function(u) {
    u <- computeSumFactors(u)
    u <- normalize(u)
    u <- runTSNE(u)
    u <- runUMAP(u)
})

labs <- c(
    sprintf("%s%% DE, 0%% DS", with(sim_pars[1:3, ], p_type * 100)),
    sprintf("%s%% DS, 0%% DE", with(sim_pars[4:6, ], unlist(map(p_dd, 3)) * 100)),
    sprintf("%s%% d(logFC), 5%% DS", with(sim_pars[7:9, ], unlist(map(map(rel_lfc, diff), 1)) * 100)))

df <- map(sim2, function(u) {
    dr <- reducedDim(u, "UMAP")
    df <- data.frame(dr, colData(u)) %>% 
        mutate_if(is.factor, as.character)
}) %>% bind_rows(.id = "i") %>% 
    mutate(id = factor(paste(cluster_id, group_id))) %>% 
    mutate_at("i", factor, labels = labs)

cols <- CATALYST:::.cluster_cols[
    seq_len(nlevels(df$id))] %>% 
    set_names(levels(df$id))

p <- ggplot(df, aes(x = X1, y = X2, col = id)) +
    facet_wrap(~ i, ncol = 3) + 
    geom_point(size = 0.4, alpha = 0.2) + 
    scale_color_manual(values = cols) +
    .prettify("void") + theme(
        panel.spacing = unit(2, "mm"),
        axis.text = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        legend.position = "none",
        strip.text = element_text(size = 8, hjust = 0))

ggsave("figures/kang/sim_pars.pdf", p,
    width = 15, height = 16, units = "cm",
    dpi = 300, useDingbats = FALSE)








