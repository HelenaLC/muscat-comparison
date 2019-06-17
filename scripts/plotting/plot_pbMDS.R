# load packages
suppressPackageStartupMessages({
    library(cowplot)
    library(dplyr)
    library(edgeR)
    library(ggplot2)
    library(muscat)
})

# source utils
source(snakemake@input$utils)

# load reference dataset
sce <- readRDS(snakemake@input$sce)

# subset to 4 clusters, 3 samples
n_cells <- table(sce$cluster_id, sce$sample_id)
n_cells <- .filter_matrix(n_cells, dim = c(nk <- 3, ns <- 4))
sce <- .filter_sce(sce, dimnames(n_cells))

# simulation parameters
de5 <- c(0.99, 0, 0.01, rep(0, 3))
defs <- list(x = sce, n_genes = nrow(sce), n_cells = 2*nk*ns*100)
sim_pars <- list(
    c(defs, list(p_type = 0.0)),
    c(defs, list(p_type = 0.01)),
    c(defs, list(p_type = 0.05)),
    
    c(defs, list(p_dd = c(0.99, 0, 0.01, rep(0, 3)))),
    c(defs, list(p_dd = c(0.95, 0, 0.05, rep(0, 3)))),
    c(defs, list(p_dd = c(0.9 , 0, 0.1,  rep(0, 3)))),
    
    c(defs, list(p_dd = de5, rel_lfc = NULL)),
    c(defs, list(p_dd = de5, rel_lfc = c(0.9,1,1.1))),
    c(defs, list(p_dd = de5, rel_lfc = c(0.8,1,1.2))))

# simulate data & aggregate
set.seed(1985)
kids <- paste0("cluster", seq_len(nk))
sims <- lapply(sim_pars, function(u) {
    sim <- do.call("simData", u)
    kids <- kids[as.numeric(u$cluster_id)]
    sim$cluster_id <- factor(kids)
    return(sim)
})
pbs <- lapply(sims, aggregateData)

# pb-level MDS plots
mains <- c(
    sprintf("%s%% type-genes, 0%% DE", c(0, 1, 5)),
    sprintf("%s%% DE, 0%% type-genes", c(1, 5, 10)),
    sprintf("%s%% \U0394(logFC), 5%% DE", c(0, 10, 20)))
ps <- lapply(seq_along(pbs), function(i) {
    p <- pbMDS(pbs[[i]]) + labs(title = mains[i], 
        x = "MDS dim. 1", y = "MDS dim. 2") +
        prettify(theme = "bw", aspect.ratio = 1) +
        theme(plot.title = element_text(size = 6))
    p$layers[[1]]$aes_params$size <- 3
    return(p)
})

# store & remove legend
lgd <- get_legend(
    ps[[1]] + theme(legend.position = "bottom") +
        guides(
            color = guide_legend(ncol = 3),
            shape = guide_legend(ncol = 2)))
ps <- lapply(ps, "+", theme(
    plot.margin = unit(c(5,rep(0,3)), "mm"), 
    legend.position = "none"))

# remove redundant y-labs
keep_ylab <- seq(1, length(ps), 3)
ps[-keep_ylab] <- lapply(ps[-keep_ylab], "+", 
    theme(axis.title.y = element_text(color = "white")))
ps <- split(ps, rep(seq_len(3), each = 3))

# round axis-limits to nearest .5 & 
# set equal within panel for comparability
n <- 1/0.5
lims <- lapply(ps, function(u) {
    xys <- lapply(u, function(v) apply(v$data[, c("MDS1", "MDS2")], 2, range))
    rng <- apply(do.call("rbind", xys), 2, range)
    rbind(floor(rng[1, ]*n)/n, ceiling(rng[2, ]*n)/n)
})
ps <- lapply(seq_along(ps), function(i)
    lapply(ps[[i]], function(p)
        p + scale_x_continuous(limits = lims[[i]][, 1]) +
            scale_y_continuous(limits = lims[[i]][, 2])))

# arrange figure
ps <- lapply(ps, function(u) plot_grid(plotlist = u, ncol = 3))
p <- plot_grid(plotlist = ps, ncol = 1,
    labels = c("a", "b", "c"), label_size = 10)
p <- plot_grid(p, lgd, ncol = 1, rel_heights = c(3, 0.2))

ggsave(snakemake@output$fig, p, 
    width = 15, height = 18, units = "cm", 
    dpi = 300, device = cairo_pdf)
