# load packages
suppressPackageStartupMessages({
    library(cowplot)
    library(ggplot2)
    library(scater)
})


source("scripts/utils.R")

# load data
sce <- readRDS("data/raw_data/kang-sim.rds")

# filter to 4 clusters & 4 samples
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, dim = c(3, 3))
sce <- .filter_sce(sce, dimnames(n_cells))

# simulation parameters
sim_pars <- list()
defs <- list(n_genes = nrow(sce), n_cells = 2*100*prod(dim(n_cells)))

sim_pars <- c(sim_pars, 
    list(c(defs, list(p_type = 0.01))),
    list(c(defs, list(p_type = 0.05))),
    list(c(defs, list(p_type = 0.10))),
    
    list(c(defs, list(p_dd = c(0.99, 0, 0.01, rep(0, 3))))),
    list(c(defs, list(p_dd = c(0.95, 0, 0.05, rep(0, 3))))),
    list(c(defs, list(p_dd = c(0.90, 0, 0.10, rep(0, 3))))),
    
    list(c(defs, list(p_dd = c(0.95, 0, 0.05, rep(0, 3)), rel_lfc = rep(1, 3)))),
    list(c(defs, list(p_dd = c(0.95, 0, 0.05, rep(0, 3)), rel_lfc = c(0.8, 1, 1.2)))),
    list(c(defs, list(p_dd = c(0.95, 0, 0.05, rep(0, 3)), rel_lfc = c(0.5, 1, 1.5)))))

# data simulation
set.seed(88)
sims <- lapply(sim_pars, function(pars) {
    sim <- do.call("simData", c(list(x = sce), pars))
    sim <- normalize(sim)
    sim <- runTSNE(sim)
    sim <- runUMAP(sim)
    sim <- runPCA(sim)
    sim <- runMDS(sim)
})

sims2 <- lapply(sims, function(u) {
    u$id <- paste0(paste0("cluster", as.numeric(u$cluster_id)), sprintf("(group%s)", u$group_id))
    return(u)
})

# plot t-SNEs
mains <- c(
    sprintf("%s%% type-genes, 0%% DE", c(1, 5, 10)),
    sprintf("%s%% DE, 0%% type-genes", c(1, 5, 10)),
    sprintf("%s%% \U0394(logFC), 5%% DE", c(0, 20, 50)))

for (dr in c("PCA", "MDS", "TSNE", "UMAP")) {
    fun <- get(paste0("plot", dr))
    labs <- paste(
        switch(dr, TSNE = "t-SNE", PCA = "PC", dr),
        switch(dr, PCA = "", "dimension"), 1:2)
    ps <- lapply(seq_along(sims), function(i) {
        fun(sims2[[i]], 
            colour_by = "id", shape_by = "group_id", 
            point_size = 0.4, point_alpha = 0.4) + 
            scale_color_manual(NULL, values = muscat:::cluster_colors) +
            scale_shape_manual(values = c(19, 19)) + guides(shape = FALSE,
                color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
            prettify(theme = "bw", aspect.ratio = 1) + 
            theme(plot.title = element_text(size = 6),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.spacing.x = unit(2, "mm")) +
            labs(title = mains[i], x = "t-SNE dimension 1", y = "t-SNE dimension 2")
    })
    lgd <- get_legend(ps[[1]] + theme(
        legend.position = "bottom",
        legend.direction = "horizontal"))
    ps <- lapply(ps, "+", theme(legend.position = "none"))
    
    keep_ylab <- seq(1, length(ps), 3)
    ps[-keep_ylab] <- lapply(ps[-keep_ylab], "+", 
        theme(axis.title.y = element_text(color = "white")))
    
    idx <- split(seq_along(ps), rep(seq_len(3), each = 3))
    ps <- lapply(idx, function(i) plot_grid(
        plotlist = ps[i], nrow = 1, align = "h", axis = "tl"))
    
    # add shared legend
    ps <- c(ps, list(lgd))
    
    p <- plot_grid(plotlist = ps, 
        ncol = 1, rel_heights = c(rep(1, 3), 0.2),
        labels = c("a", "b", "c"), label_size = 10)
    
    # save figure to .pdf
    ggsave(sprintf("figures/sim_pars_%s.pdf", dr), p,
        width = 15, height = 18, units = "cm",
        dpi = 300, device = cairo_pdf)
}

# t-SNE only
ps <- lapply(seq_along(sims), function(i) {
    plotTSNE(sims2[[i]], 
        colour_by = "id", shape_by = "group_id", 
        point_size = 0.2, point_alpha = 0.4) + 
        scale_color_manual(NULL, values = muscat:::cluster_colors) +
        scale_shape_manual(values = c(19, 19)) + guides(shape = FALSE,
            color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
        prettify(theme = "bw", aspect.ratio = 1) + 
        theme(plot.title = element_text(size = 6),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(rep(0, 4), "cm"),
            axis.ticks.length = unit(0, "pt"))
})
p <- plot_grid(plotlist = ps, ncol = 3, align = "h", axis = "tl")
ggsave("/users/helena/desktop/fig1c.pdf", p,
    width = 10, height = 10, units = "cm",
    dpi = 300, useDingbats = FALSE)

