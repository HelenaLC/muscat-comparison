suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(magrittr)
    library(purrr)
    library(reshape2)
})

print(args)
# load performance plots
ps <- lapply(args$ggp, readRDS)
lgd1 <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()))

# load simulation parameters
sim_pars <- lapply(args$sim_pars, yaml::read_yaml)
n <- sapply(c("nk", "ns", "nc"), function(u) unlist(map(sim_pars, u)))
n <- with(as.data.frame(n), nc / (2 * nk * ns))
ss <- map(map(sim_pars, "probs"), 2) 
ss <- sapply(seq_along(ss), function(i) n[i] * ss[[i]]) 
colnames(ss) <- gsub(".json", "", basename(args$sim_pars))

gg_df <- melt(ss) %>% 
    set_colnames(c("sample_id", "sim_id", "n")) %>% 
    replicate(n = 2, simplify = FALSE) %>% 
    bind_rows(.id = "group_id") %>% 
    mutate_at("group_id", function(u) factor(u, labels = c("A", "B"))) %>% 
    dplyr::mutate(sample_id = factor(paste0(group_id, sample_id)))

p0 <- ggplot(gg_df, aes(x = sample_id, y = n, fill = sample_id)) +
    facet_wrap("sim_id", nrow = 1) +
    geom_bar(size = 1, width = 0.5, stat = "identity") +
    scale_fill_manual(values = c("darkblue", "royalblue", "lightblue", "tomato", "orange", "gold")) +
    scale_y_continuous("expected\nnb. of cells", limits = c(0, 125)) +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE,
        override.aes = list(size = 0))) +
    .prettify("bw") + theme(
        aspect.ratio = NULL,
        legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank())

lgd2 <- get_legend(p0)
p0 <- p0 + theme(legend.position = "none")

# arrange plots
p <- plot_grid(
    plotlist = c(list(p0), ps),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 2, 2),
    labels = c("", "a", "b"),
    label_size = 10,
    label_fontface = "bold")
l <- plot_grid(lgd1, lgd2, rel_widths = c(2.4, 1))
p <- plot_grid(p, l, ncol = 1, 
    align = "v", axis = "r",
    rel_heights = c(8, 1))

ggsave(args$fig, p,
    width = 15, height = 12.5, units = "cm",
    dpi = 300, useDingbats = FALSE)





