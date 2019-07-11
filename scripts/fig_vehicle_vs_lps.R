source(file.path("scripts", "utils.R"))
dir <- "/users/helena/dropbox/phd/manuscript/fig4"

# load packages
suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(magrittr)
    library(muscat)
    library(SingleCellExperiment)
    library(purrr)
})

# load data & results
sce <- readRDS(file.path("MAGL", "output", "MAGL-SCE.rds"))
#sce <- readRDS("/Volumes/Shared/data/seq/calini_scrnaseq_output/MAGL_mouse/output/MAGL-SCE.rds")

# construct data.frame for plotting
ids <- c("cluster_id", "group_id", "sample_id")
df <- data.frame(colData(sce), do.call(cbind, reducedDims(sce)))

# store IDs & nb. of levels
nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

# specify color palettes
cluster_id_pal <- CATALYST:::.cluster_cols[seq_len(nk)]
sample_id_pal <- CATALYST:::.cluster_cols[seq_len(ns) + nk]
group_id_pal <- c("royalblue", "orange")

# color legends ----------------------------------------------------------------
for (id in ids) {
    fn <- sprintf(file.path(dir, "lgd-%s.pdf"), id)
    pal <- get(paste0(id, "_pal"))
    foo <- ggplot(df, aes_string("ident", "ident", col = id)) + 
        scale_color_manual(values = pal) + geom_point() + .prettify()
    ggsave(fn, get_legend(foo), 
        width = 3, height = 4, units = "cm", 
        dpi = 300, useDingbats = FALSE)
}

# relative cluster abundances --------------------------------------------------
p <- table(sce$cluster_id, sce$sample_id) %>% 
    prop.table(2) %>% unclass %>% melt %>% 
    set_colnames(c("cluster_id", "sample_id", "frequency")) %>% 
    dplyr::mutate(group_id = sce$group_id[match(.$sample_id, sce$sample_id)]) %>% 
    ggplot(aes(x = sample_id, y = frequency, fill = cluster_id)) +
    geom_bar(stat = "identity", width = 1, col = "white", show.legend = FALSE) + 
    facet_wrap("group_id", scales = "free_x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = cluster_id_pal) +
    .prettify() + theme(aspect.ratio = 2, panel.spacing = unit(2, "mm"), 
        strip.background = element_blank(), strip.text = element_blank(),
        panel.grid = element_blank(), axis.line = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file.path(dir, "cluster_freqs.pdf"), p,
    width = 12, height = 12, units = "cm",
    dpi = 300, useDingbats = FALSE)   

# UMAP colored by cluster, sample & group ID -----------------------------------
.plot_dr <- function(df, col, pal)
    ggplot(df, aes_string(x = "UMAP_1", y = "UMAP_2", col = col)) +
    geom_point_rast(size = 0.1, alpha = 0.1, show.legend = FALSE) + 
    scale_color_manual(values = pal) + 
    theme_void() + theme(aspect.ratio = 1)

for (id in ids) {
    p <- .plot_dr(df, id, get(paste0(id, "_pal")))
    fn <- sprintf(file.path(dir, "umap-%s.pdf"), id)
    ggsave(fn, p, 
        width = 10, height = 10, units = "cm", 
        dpi = 300, useDingbats = FALSE)
}

# pseudobulk-level MDS plot ----------------------------------------------------
pb <- aggregateData(sce)
mds <- pbMDS(pb)


