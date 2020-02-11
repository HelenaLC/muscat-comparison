# load config & source utils
config <- yaml::read_yaml("config.yaml")
source(config$utils)

# load packages
suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(reshape2)
})

# load method metadata
md <- read.csv(config$mids, row.names = 2)[, -1, drop = FALSE]
md$id <- rownames(md)

lvls <- factor(seq_len(3), labels = c("good", "intermediate", "poor"))
error_control <- TPR <- complex_design <- speed <- setNames(rep("good", nrow(md)), md$id)
logFC_est <- setNames(rep(NA, nrow(md)), md$id)

# error control
error_control[c("edgeR.sum(scalecpm)", "MM-vst", "scDD.logcounts")] <- "intermediate"
error_control[c("scDD.vstresiduals", md$id[md$type == "ad"])] <- "poor"

# TPR
TPR[c("MM-vst", "scDD.logcounts", "edgeR.sum(scalecpm)", "MM-nbinom")] <- "intermediate"

# logFC estimation
logFC_est[c(grep("edgeR", md$id, value = TRUE), "limma-voom.sum(counts)")] <- "good"
logFC_est[c(grep("trend", md$id, value = TRUE), "MM-dream")] <- "intermediate"

# design
complex_design[md$id[md$type != "pb"]] <- "poor"
complex_design[c("MAST.logcounts", grep("MM", md$id, value = TRUE))] <- "intermediate"

# runtime
speed[grep("MM", md$id, value = TRUE)] <- "poor"
speed[md$id[md$type %in% c("ad", "scdd")]] <- "intermediate"

# plot summary heatmap --------------------------------------------------------
df <- data.frame(md, error_control, TPR, logFC_est, complex_design, speed) %>% 
    melt(id.vars = "id") %>% mutate_at("variable", factor) %>% 
    mutate_at("value", factor, levels = c(c("pb", "mm", "ad", "mast", "scdd"), levels(lvls)))

pal <- c(.typ_cols, setNames(c("skyblue", "gold", "tomato"), lvls))

ys <- rev(levels(df$variable))
xs <- dplyr::filter(df, variable != "type") %>% group_by(id) %>% 
    dplyr::mutate(score = mean(match(value, rev(lvls), nomatch = 0))) %>% 
    dplyr::slice(1) %>% ungroup %>% 
    dplyr::arrange(desc(score)) %>% pull(id)

p <- ggplot(df, aes(x = id, y = variable, fill = value)) +
    geom_point() +
    geom_tile(data = df[df$variable == "type", ], height = 0.4, size = 0.5, col = "white") +
    geom_tile(data = df[df$variable != "type", ], size = 0.5, col = "white") +
    scale_fill_manual(NULL, na.value = "lightgrey", values = pal,
        labels = c(.typ_labs, "good", "inter-\nmediate", "poor")) +
    guides(fill = guide_legend(override.aes = list(size = 0.1), 
        byrow = FALSE, nrow = nlevels(md$type))) +
    coord_cartesian(clip = "off") +
    scale_x_discrete(limits = xs, expand = c(0, 0)) +
    scale_y_discrete(limits = ys, expand = c(0, 0)) +
    .prettify("void") + theme(
        aspect.ratio = nlevels(df$variable) / nrow(md),
        plot.margin = unit(rep(0, 4), "mm"),
        legend.justification = "top",
        legend.margin = margin(-1.3,0,0,5, "mm"),
        legend.key.height = unit(4.5, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(1, "mm"), 
        axis.text.x = element_text(size = 6, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6, hjust = 1))

ggsave(file.path("figures", "summary.pdf"), p,
    width = 15, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)
