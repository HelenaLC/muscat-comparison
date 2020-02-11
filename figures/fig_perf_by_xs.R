suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(magrittr)
    library(purrr)
    library(RColorBrewer)
    library(reshape2)
})

for (x in c("s", "g")) {

# load performance plots
fns <- sprintf("%s-perf_by_%ss.rds", config$dids, x)
ps <- lapply(file.path("plots", fns), readRDS)
lgd1 <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()))

# load simulation parameters
sids <- unique(unlist(map(map(ps, "data"), "sid")))
fns <- paste0(config$sim_pars, sids, ".json")
sim_pars <- lapply(fns, yaml::read_yaml)
names(ns) <- ns <- c("nk", "ns", "nc")
ns <- lapply(ns, function(u) unlist(map(sim_pars, u)))
nc <- with(as.data.frame(ns), nc / (2 * nk * ns))
ss <- map(map(sim_pars, "probs"), ifelse(x == "s", 2, 3)) 
ss <- sapply(seq_along(ss), function(i) nc[i] * ss[[i]])
colnames(ss) <- gsub(".json", "", basename(fns))

if (x == "s") {
    gg_df <- melt(ss) %>% 
        set_colnames(c("sample_id", "sim_id", "n")) %>% 
        replicate(n = 2, simplify = FALSE) %>% 
        bind_rows(.id = "group_id")
} else {
    gg_df <- melt(ss) %>% 
        set_colnames(c("group_id", "sim_id", "n")) %>% 
        replicate(n = ns$ns[1], simplify = FALSE) %>% 
        bind_rows(.id = "sample_id")
}
gg_df <-  mutate_at(gg_df, "group_id", 
    function(u) factor(u, labels = c("A", "B"))) %>% 
    dplyr::mutate(sample_id = factor(paste0(group_id, sample_id)))

ns <- nlevels(gg_df$sample_id)/2
cols <- lapply(c("Blues", "Reds"), function(pal)
    brewer.pal(9, pal)[seq(7, 9-2*ns, -2)])

p0 <- ggplot(gg_df, aes(x = sample_id, y = n, fill = sample_id)) +
    facet_wrap("sim_id", nrow = 1) +
    geom_bar(size = 1, width = 0.5, stat = "identity") +
    scale_fill_manual(values = unlist(cols)) +
    scale_y_continuous("expected\nnb. of cells", expand = c(0, 0), 
        limits = c(0, ceiling(max(gg_df$n)/50)*50)) +
    guides(fill = guide_legend(ncol = ns, byrow = TRUE,
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

fn <- file.path("figures", sprintf("perf_by_%ss.pdf", x))
ggsave(fn, p, 
    width = 15, height = 12.5, units = "cm",
    dpi = 300, useDingbats = FALSE)

}