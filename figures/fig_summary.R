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
md <- read.csv(config$mids, row.names = 1)
rownames(md) <- mids <- md$id

scores <- rep("good", nrow(md))
values <- c("TPR", "FDR", "error\ncontrol", "logFC\nestimation", "complex\ndesign", "speed")
df <- replicate(length(values), scores, simplify = FALSE)
df <- data.frame(df, stringsAsFactors = FALSE)
rownames(df) <- mids; colnames(df) <- values

df[grepl("MM-vst|dream|scDD.logcounts|scalecpm", mids), "TPR"] <- "intermediate"

df[grepl("AD|scDD", mids), "FDR"] <- "poor"
df[grepl("AD-gid.log|MAST|MM-vst", mids), "FDR"] <- "intermediate"

df[grepl("AD", mids), "error\ncontrol"] <- "poor"
df[grepl("MM-(d|v)", mids), "error\ncontrol"] <- "intermediate"

df[grepl("AD|MAST|MM", mids), "speed"] <- "intermediate"
df[grepl("MM-nbinom", mids), "speed"] <- "poor"

df[grepl("AD|scDD", mids), "complex\ndesign"] <- "poor"

df[grepl("log|vst|MM", mids), "logFC\nestimation"] <- "intermediate"
df[grepl("AD|scDD|MAST", mids), "logFC\nestimation"] <- "NA"

gg_df <- cbind(md, df) %>% 
    melt(id.vars = c("id", "type")) %>% 
    mutate_at("value", factor, 
        levels = c("good", "intermediate", "poor", "NA"),
        labels = c("good", "inter-\nmediate", "poor", "NA")) %>% 
    mutate_at("type", factor, levels = names(.typ_labs), labels = .typ_labs)

ys <- c(rev(values), "method\nclass")
xs <- group_by(gg_df, id) %>% 
    mutate(score = mean(match(value, setdiff(rev(levels(value)), "NA"), nomatch = 0))) %>% 
    dplyr::slice(1) %>% ungroup %>% arrange(desc(score)) %>% pull(id)

pal <- c(.typ_cols, c("skyblue", "gold", "tomato"), "lightgrey")
names(pal) <- c(.typ_labs, levels(gg_df$value))

p <- ggplot(gg_df, aes(x = id, y = variable, fill = value)) + 
    geom_tile(height = 1, size = 0.5, col = "white") +
    geom_tile(aes(fill = type, y = length(values) + 1), 
        height = 0.5, size = 0.5, col = "white") +
    scale_fill_manual(NULL, values = pal, 
        limits = c("AD", "MAST", "scDD", "PB", "MM", levels(gg_df$value))) +
    guides(fill = guide_legend(override.aes = list(size = 0.2),
        byrow = FALSE, nrow = nlevels(gg_df$type))) +
    scale_x_discrete(limits = xs, expand = c(0, 0)) +
    scale_y_discrete(limits = ys, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    .prettify("void") + theme(
        plot.margin = unit(rep(0, 4), "mm"),
        aspect.ratio = (length(values)+0.5)/length(mids),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(1, "mm"), 
        axis.text.x = element_text(size = 6, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6, hjust = 1),
        legend.justification = "top",
        legend.margin = margin(0,0,0,0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.height = unit(6, "mm"),
        legend.key.width = unit(3, "mm"))

ggsave(file.path("figures", "summary_heatmap.pdf"), p,
    width = 15, height = 6.5, units = "cm",
    dpi = 300, useDingbats = FALSE)
