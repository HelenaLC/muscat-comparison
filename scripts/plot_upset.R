suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(purrr)
    library(UpSetR)
})

#args <- list(res = list.files("results", "kang,d[a-z][0-9]+,", full.names = TRUE))
res <- mutate(.read_res(args$res), 
    hit = paste(gene, cluster_id, sid, i, sep = ";"))

n_dd <- res %>% 
    filter(mid == .$mid[1]) %>% 
    group_by(sid, i) %>% 
    summarize(n_dd = sum(is_de)) %>% 
    acast(sid ~ i, value.var = "n_dd")

# get gene-cluster combinations by method at FDR 5%
top <- group_by(res, sid, i) %>% do(
    dplyr::filter(., p_adj.loc < 0.05) %>% 
        dplyr::arrange(p_val) %>% 
        group_by(mid, add = TRUE) %>%
        dplyr::slice(seq_len(n_dd[.$sid[1], .$i[1]])) %>% 
        summarize(hit = list(hit))) %>% 
    group_by(mid) %>% 
    summarize(hit = list(purrr::reduce(hit, c))) %>% 
    mutate_at("mid", as.character)

# get intersections & gene metadata
df_bars <- fromList(set_names(top$hit, top$mid)) %>% 
    dplyr::mutate(
        code = apply(.[top$mid], 1, paste, collapse = ""),
        degree = apply(.[top$mid], 1, sum),
        hit = unique(unlist(top$hit))
    ) %>% {
        m <- match(.$hit, res$hit)
        mutate(., sid = res$sid[m], i = res$i[m], cat = res$category[m])
    } %>% add_count(code) %>% group_by(code) %>% ungroup

# order by degree & select most frequent interactions
m <- match(unique(df_bars$code), df_bars$code)
keep <- pull(top_n(df_bars[m, ], 40, n), "code")
df_bars <- dplyr::filter(df_bars, code %in% keep)
m <- match(unique(df_bars$code), df_bars$code)
o <- order(df_bars$degree[m], -df_bars$n[m])

# construct data.frame of unique method intersections
df_grid <- melt(df_bars, variable.name = "mid",
    id.var = setdiff(names(df_bars), top$mid)) %>% 
    group_by(code, mid) %>% slice(1)

# get method type annotation
df_anno <- read.csv(config$mids)
df_anno$type <- factor(df_anno$type, 
    levels = names(.typ_cols))

# plotting ---------------------------------------------------------------------

df_bars_ee <- filter(df_bars, cat == "ee")
df_bars_dx <- filter(df_bars, cat != "ee")

x_lims <- df_bars$code[m][o]
y_lims <- rev(levels(res$mid))

y_max_ee <- ceiling(max(df_bars_ee$n)/1e3)*1e3
y_max_dx <- y_max_ee / 10

# shared aesthetics
thm <- .prettify("classic") + theme(
    aspect.ratio = NULL,
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(0,0,0,0, "mm"))

# barplot of intersection sizes
bars_ee <- ggplot(df_bars_ee, aes(x = code)) +
    stat_count(aes(fill = cat)) + 
    scale_fill_manual("gene\ncategory", limits = levels(df_bars$cat),
        values = .cat_cols, labels = toupper(levels(df_bars$cat))) +
    guides(fill = guide_legend(
        override.aes = list(col = "white", size = 0.1))) +
    scale_x_discrete(limits = x_lims, expand = c(0,0)) +
    scale_y_continuous(limits = c(0, y_max_ee), expand = c(0,0)) + 
    thm + theme(
        plot.margin = unit(c(1,0,1,1), "mm"),
        panel.grid.major.x = element_blank())
bars_dx <- ggplot(df_bars_dx, aes(x = code)) +
    stat_count(aes(fill = cat), show.legend = FALSE) + 
    scale_x_discrete(limits = x_lims, expand = c(0,0)) +
    scale_y_continuous(limits = c(0, y_max_dx), expand = c(0,0)) + 
    scale_fill_manual(NULL, values = .cat_cols, 
        labels = toupper(names(.cat_cols))) + 
    thm + theme(
        plot.margin = unit(c(1,0,1,1), "mm"),
        panel.grid.major.x = element_blank())

# method intersection grid
grid <- ggplot(df_grid, aes(x = code, y = mid, color = factor(value))) +
    scale_color_manual(values = c("0" = "grey90", "1" = "black")) +
    scale_x_discrete(limits = x_lims) +
    scale_y_discrete(limits = y_lims) +
    geom_point(shape = 16, size = 1) +
    geom_path(aes(group = code), size = 0.2,
        data = filter(df_grid, value != 0)) + 
    annotate("rect", alpha = 0.08, xmin = 0.5, xmax = Inf, 
        ymin = seq(0.5,nrow(top),2), ymax = seq(1.5, nrow(top)+1, 2)) +
    thm + theme(
        legend.position = "none", 
        panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")) 

# method class annotation
anno <- ggplot(df_anno, aes(x = 0, y = id, fill = type)) +
    scale_fill_manual("method\nclass", 
        labels = .typ_labs, values = .typ_cols) + 
    scale_y_discrete(limits = y_lims) +
    geom_tile(col = "white", size = 0.1) + coord_fixed(1) + 
    thm + theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"))

# get legends of DD category & method annotations
lgd_bars <- get_legend(bars_ee)
lgd_anno <- get_legend(anno + theme(legend.justification = "top")
bars_ee <- bars_ee + theme(legend.position = "none")
anno <- anno + theme(legend.position = "none")

# arrange plots
foo <- ggplot() + theme_nothing()
p <- plot_grid(
    plot_grid(bars_ee, bars_dx, grid, 
        ncol = 1, align = "v", axis = "lr"),
    plot_grid(foo, foo, anno, ncol = 1),
    plot_grid(lgd_bars, lgd_anno,
        ncol = 1, rel_heights = c(2, 1),
        align = "v", axis = "t"),
    nrow = 1, rel_widths = c(15, 0.5, 2))

ggsave(args$fig, p, dpi = 300,
    width = 15, height = 9, units = "cm")

