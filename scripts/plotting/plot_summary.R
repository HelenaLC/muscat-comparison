# load packages
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(dplyr)
    library(ggplot2)
    library(grid)
    library(reshape2)
})

# source utils
source("scripts/utils.R")

lo <- 0; md <- 1; hi <- 2
cols <- c("tomato", "gold", "lightskyblue")

# type 1 error control
res <- readRDS("results/kang/nill_df.rds") %>% 
    group_by(replicate, method) %>% 
    summarise(ks.test(p_val, "punif", 0, 1)$p.value) %>% 
    acast(replicate ~ method) %>% 
    apply(2, "<", 0.05) %>% 
    colMeans

df <- read.csv("metadata/method_ids.csv", 
    row.names = 1, stringsAsFactors = FALSE)

df$error_control <- md
df[names(res)[res == 0], "error_control"] <- lo
df[names(res)[res == 1], "error_control"] <- hi

df$TPR <- hi
df[df$type == "mast", "TPR"] <- md
df[df$type %in% c("ad", "scdd"), "TPR"] <- lo

df$logFC_est <- hi
df[df$type == "mm", "logFC_est"] <- md
df[df$id %in% c("limma-trend.mean(normcounts)", "limma-trend.mean(logcounts)"), "logFC_est"] <- md
df[df$id == "limma-trend.sum(logcounts)", "logFC_est"] <- lo

df$complex_design <- hi
df[df$type %in% c("ad", "mast", "scdd"), "complex_design"] <- lo

df$speed <- hi
df[df$type == "mm", "speed"] <- lo
df[!df$type %in% c("pb", "mm"), "speed"] <- md

o <- df %>% 
    select(-c("id", "type")) %>% 
    rowMeans %>% 
    sort(decreasing = TRUE) %>% 
    names

gg_df <- melt(select(df, -"type"), id.var = "id")
gg_df <- melt(df, id.var = c("id", "type"))
p <- ggplot(gg_df, aes(x = id, y = variable, fill = factor(value))) +
    geom_tile(col = "white", width = 1, height = 1, size = 0.5) +
    scale_fill_manual(NULL, values = cols, breaks = c(hi, md, lo),
        labels = c("good", "inter-\nmediate", "poor")) +
    guides(fill = guide_legend(override.aes = list(size = 0))) +
    prettify(theme = "void") +
    scale_x_discrete(limits = o) +
    scale_y_discrete(limits = rev(levels(gg_df$variable))) +
    annotate("rect", xmin = 0.5, xmax = 14.5, ymin = 0.5, ymax = nlevels(gg_df$variable) + 0.5,
        fill = NA, col = "black", lty = 2, size = 0.25) +
    coord_fixed() + theme(
        legend.margin = margin(0,1,0,0, "mm"),
        legend.key.size = unit(3, "mm"),
        legend.justification = "right",
        legend.position = "bottom",
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length=unit(1, "mm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            #face = "bold", color = method_colors[o]),
        axis.text.y = element_text(hjust = 1))

ggsave("figures/method_summary.pdf", p,
    width = 15, height = 7, units = "cm",
    dpi = 300, useDingbats = FALSE)

