# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(cowplot)
        library(dplyr)
        library(muscat)
        library(purrr)
        library(scater)
    }))

# source utils
source(snakemake@input$utils)
set.seed(1)

# 1 group, x% type, no DE ------------------------------------------------------
sce <- readRDS("/users/helena/dropbox/phd/projects/muscat-comparison/data/raw_data/kang-sim.rds")
sce <- scater::normalize(sce)
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, dim = c(3, 3))
sce <- .filter_sce(sce, dimnames(n_cells))

p_type <- c(0.01, 0.05, 0.1)
names(p_type) <- p_type
dr0 <- .run_dr(sce)
dr <- lapply(p_type, function(p) {
    simData(sce, 
        n_genes = nrow(sce), n_cells = 16*200, p_dd = diag(6)[1, ], 
        probs = list(NULL, NULL, c(0, 1)), p_type = p) %>% 
        normalize %>% .run_dr
}) %>% c(ref = dr0)
p <- .plot_dr(dr, color_by = "cluster_id", lab = "type-genes")
ggsave("figures/tsne_ptype.pdf", p,
    width = 15, height = 4.5, unit = "cm",
    dpi = 300, useDingbats = FALSE)

# 2 groups, no type, x% DE -----------------------------------------------------
sce <- readRDS("/users/helena/dropbox/phd/projects/muscat-comparison/data/raw_data/kang-sce.rds")
sce <- scater::normalize(sce)
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, n = 200)
sce <- .filter_sce(sce, dimnames(n_cells))
sce <- normalize(sce)

p_de <- c(0.01, 0.05, 0.1)
names(p_de) <- p_de
dr0 <- .run_dr(sce)
dr <- lapply(p_de, function(p) {
    simData(sce, 
        n_genes = nrow(sce), n_cells = 16*200, 
        p_dd = c(1-p,0,p,rep(0,3)), p_type = 0) %>% 
        normalize %>% .run_dr
}) %>% c(ref = dr0)
p1 <- .plot_dr(dr, color_by = "cluster_id", lab = "DE")
p2 <- .plot_dr(dr, color_by = "group_id", lab = "DE")
p <- plot_grid(p1, p2, nrow = 2,
    axis = "tl", align = "hv")
ggsave("figures/tsne_pde.pdf", p,
    width = 15, height = 8.5, unit = "cm",
    dpi = 300, useDingbats = FALSE)

# 1 group, x% type, no DE ------------------------------------------------------
sce <- readRDS("/users/helena/dropbox/phd/projects/muscat-comparison/data/raw_data/kang-sim.rds")
sce <- scater::normalize(sce)
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, dim = c(3, 3))
sce <- .filter_sce(sce, dimnames(n_cells))

rel_lfc <- list(rep(1, 3), c(1.2,1,0.8), c(2,1,0.5))
names(rel_lfc) <- c(0, 0.2, 0.5)
dr0 <- .run_dr(sce)
dr <- lapply(rel_lfc, function(x) {
    simData(sce, 
        n_genes = nrow(sce), n_cells = 16*200, 
        p_dd = c(.9,0,.1,0,0,0), lfc = 3, rel_lfc = x) %>% 
        normalize %>% .run_dr
}) %>% c(ref = dr0)
p <- .plot_dr(dr, color_by = "cluster_id", lab = "rel. diff. logFC")
ggsave("figures/tsne_rel_lfc.pdf", p,
    width = 15, height = 4.5, unit = "cm",
    dpi = 300, useDingbats = FALSE)
