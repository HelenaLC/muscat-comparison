source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(purrr)
})

#fns <- list.files("~/documents/kang", "d[a-z][0-9]+;", full.names = TRUE)
df <- lapply(snakemake@input$res, readRDS) %>% 
    map("tbl") %>% 
    map(mutate_if, is.factor, as.character) %>% 
    bind_rows %>% 
    dplyr::filter(!is.na(sim_lfc) & !is.na(est_lfc)) %>% 
    rename(method = mid) %>% 
    mutate_at("method", function(u) droplevels(
        factor(u, levels = names(.meth_cols)))) %>% 
    mutate_at("sid", factor, 
        levels = paste0(c("ds", "dp", "dm", "db"), "10"),
        labels = c("DS", "DP", "DM", "DB"))

# dowsample to 2k cluster-gene combinations
# per method, simulation & relicate
set.seed(29)
sub <- ungroup(sample_n(group_by(df, method, sid, i), 2e3))

# cut 0.1% quantiles
rng <- quantile(c(sub$sim_lfc, sub$est_lfc), c(1e-3, 1-1e-3))
sub <- mutate_at(sub, c("sim_lfc", "est_lfc"), function(u) {
    u[u < rng[1]] <- rng[1]; u[u > rng[2]] <- rng[2]; return(u) })

p <- ggplot(sub, aes(x = sim_lfc, y = est_lfc, col = as.logical(is_de))) +
    facet_grid(vars(sid), vars(method)) +
    geom_point_rast(size = 1.2, alpha = 0.4, raster.dpi = 100) +
    scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "tomato")) +
    guides(color = guide_legend("differential", override.aes = list(size = 3, alpha = 1))) +
    scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,3), expand = c(0, 1)) +
    scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3), expand = c(0, 1)) +
    labs(x = "simulated logFC", y = "estimated logFC") +
    .prettify("bw") + theme(
        strip.text = element_text(size = 4),
        panel.spacing = unit(1, "mm"),
        legend.position = "bottom")

ggsave(snakemake@output$fig, p,
    width = 15, height = 11, units = "cm",
    dpi = 300, useDingbats = FALSE)
