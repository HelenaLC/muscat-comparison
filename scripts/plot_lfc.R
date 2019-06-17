source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(purrr)
})

df <- lapply(snakemake@input$res, readRDS) %>% 
    map("tbl") %>% bind_rows %>% 
    filter(!is.na(est_lfc)) %>% 
    mutate_at("mid", factor, levels = names(.meth_cols)) %>% 
    mutate_at("mid", droplevels) %>% 
    rename(method = mid)

ps <- group_by(df, sid) %>% {set_names(group_split(.), group_keys(.)[[1]])} %>% lapply(function(u) 
    ggplot(u, aes(x = sim_lfc, y = est_lfc, col = as.logical(is_de))) +
        facet_wrap(~ method, ncol = 5) +
        geom_point(size = 0.2, alpha = 0.2) +
        scale_color_manual(values = c("FALSE" = "royalblue", "TRUE" = "tomato")) +
        guides(color = guide_legend("DD", override.aes = list(size = 3, alpha = 1))) +
        scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, 3), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-7, 7), breaks = seq(-6, 6, 3), expand = c(0, 0)) +
        labs(x = "simulated logFC", y = "estimated logFC") +
        .prettify("bw") + theme(
            strip.text = element_text(size = 5),
            panel.spacing = unit(0.1, "cm"),
            legend.position = "bottom"))

lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))
names(ps) <- toupper(gsub("[0-9]", "", names(ps)))
ps <- ps[c("DS", "DP", "DM", "DB")]

p <- plot_grid(
    plotlist = c(ps, list(lgd)), 
    ncol = 1,
    rel_heights = c(rep(1, 4), 0.2),
    labels = names(ps),
    label_size = 10,
    label_fontface = "bold")

ggsave(snakemake@output$fig, p,
    width = 15, height = 28, units = "cm",
    dpi = 300, useDingbats = FALSE)
