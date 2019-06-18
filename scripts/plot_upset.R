source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(purrr)
})

#fns <- list.files("results/kang", "d[a-z]10;", full.names = TRUE)
res <- lapply(snakemake@input$res, readRDS) %>%
    map("tbl") %>% bind_rows %>% 
    mutate(hit = paste(gene, cluster_id, sid, i, sep = ";")) %>% 
    mutate_at("category", plyr::revalue, c("de" = "ds"))

n_dd <- res %>% 
    dplyr::filter(mid == .$mid[1]) %>% 
    group_by(sid, i) %>% 
    summarize(n_dd = sum(is_de)) %>% 
    acast(sid ~ i, value.var = "n_dd")

top <- res %>% group_by(sid, i) %>% do({
    n <- n_dd[.$sid[1], .$i[1]]
    arrange(., p_val) %>%
        group_by(mid, add = TRUE) %>%
        slice(seq_len(n)) %>% 
        summarize(hit = list(hit))
}) %>% group_by(mid) %>% summarize(hit = list(reduce(hit, c)))

.cat_cols <- c("royalblue", "cornflowerblue", "red3", "tomato", "orange", "gold")
names(.cat_cols) <- levels(res$category)

df <- UpSetR::fromList(set_names(top$hit, top$mid)) %>% mutate(
    code = apply(.[top$mid], 1, paste, collapse = ""),
    degree = apply(.[top$mid], 1, sum),
    hit = unique(unlist(top$hit))) %>% {
        m <- match(.$hit, res$hit)
        mutate(., sid = res$sid[m], i = res$i[m], cat = res$category[m])
    } %>% add_count(code) %>% group_by(code) %>% 
    mutate(p_true = 100 * sum(cat %in% c("ds", "dp", "dm", "db")) / n) %>% 
    ungroup

m <- match(unique(df$code), df$code)
keep <- pull(top_n(df[m, ], 50, n), "code")
df <- filter(df, code %in% keep)
m <- match(unique(df$code), df$code)
o <- order(df$degree[m], -df$n[m])

#/ (sum(n_dd) * nrow(top))
max <- ceiling(max(df$n)/500)*500
p1 <- ggplot(df, aes(x = code)) +
    coord_cartesian(clip = 'off') +
    geom_bar(aes(y = ..count.., fill = cat)) + 
    geom_point(shape = 17, size = 1, col = "lightgrey", 
        aes(y = p_true * max / 100)) +
    scale_x_discrete(limits = df$code[m][o]) +
    scale_y_continuous(limits = c(0, max), expand = c(0,0),
        sec.axis = sec_axis(~./max(.), breaks = seq(0, 1, 0.2))) +
    scale_fill_manual(NULL, values = .cat_cols, 
        labels = function(u) toupper(u)) +
    .prettify("classic") + theme(
        legend.margin = margin(0,0,0,0, "mm"),
        plot.margin = unit(c(2,2,2,0), "mm"),
        legend.key.size = unit(4, "mm"),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = NULL)

lgd <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

dfm <- melt(df, id.var = setdiff(names(df), top$mid), variable.name = "method")
p2 <- ggplot(dfm, aes(x = code, y = method, color = factor(value))) +
    scale_x_discrete(limits = df$code[m][o]) +
    scale_y_discrete(limits = rev(top$mid)) +
    scale_color_manual(values = c("0" = "grey92", "1" = "black")) +
    geom_point(shape = 16, size = 1) +
    geom_path(size = 0.2, data = filter(dfm, value != 0), aes(group = code)) +
    annotate("rect", alpha = 0.08, xmin = 0.5, xmax = Inf, 
        ymin = seq(0.5,nrow(top),2), ymax = seq(1.5,nrow(top)+1,2)) +
    .prettify("classic") + theme(
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = NULL) 

anno <- read.csv(snakemake@config$mids)
anno$type <- factor(anno$type, levels = unique(anno$type[match(top$mid, anno$id)]))
.typ_cols <- setNames(RColorBrewer::brewer.pal(nlevels(anno$type), "Set2"), levels(anno$type))

p3 <- ggplot(anno, aes(x = 0, y = id, fill = type)) +
    scale_fill_manual("type", values = .typ_cols,
        labels = c(ad = "AD", mast = "MAST", pb = "PB", scdd = "scDD")) +
    geom_tile(col = "white") + 
    scale_y_discrete(limits = rev(top$mid)) +
    .prettify("classic") + theme(aspect.ratio = NULL,
        legend.margin = margin(0,0,0,0, "mm"),
        plot.margin = unit(c(0,2,0,0), "mm"),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

p <- plot_grid(
    plot_grid(p1, p2, ncol = 1, 
        align = "v", axis = "lr", 
        rel_heights = c(3, 2)),
    plot_grid(lgd, p3, ncol = 1, 
        align = "v", axis = "l", 
        rel_heights = c(3, 2)),
    nrow = 1, rel_widths = c(6, 1))

ggsave(snakemake@output$fig, p, dpi = 300,
    width = 15, height = 10, units = "cm")
