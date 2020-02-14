suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(purrr)
})

#fns <- list.files("~/documents/kang", "d[a-z][0-9]+;", full.names = TRUE)
res <- .read_res(args$res, wcs$inc)
mids <- levels(res$mid)

df <- rename(res, method = mid) %>%
    dplyr::filter(!(is.na(sim_lfc) | is.na(est_lfc))) %>% 
    mutate_at("method", droplevels) %>% 
    mutate_at("sid", factor, 
        levels = paste0(c("de", "dp", "dm", "db"), "10"),
        labels = c("DE", "DP", "DM", "DB"))

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
        strip.text = element_text(size = 3),
        panel.spacing = unit(1, "mm"),
        legend.position = "bottom")

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    width = 15, height = 10, units = "cm",
    dpi = 300, useDingbats = FALSE)
