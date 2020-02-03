suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(purrr)
    library(UpSetR)
})

#args <- list(res = list.files("~/projects/portmac/results", "kang,d[a-z][0-9]+,", full.names = TRUE))
res <- .read_res(args$res) %>% 
    dplyr::mutate(hit = paste(gene, cluster_id, sid, i, sep = ";"))

n_dd <- res %>% 
    dplyr::filter(mid == .$mid[1]) %>% 
    group_by(sid, i) %>% 
    summarize(n_dd = sum(is_de)) %>% 
    acast(sid ~ i, value.var = "n_dd")

top <- group_by(res, sid, i) %>% do(
    dplyr::filter(., p_adj.loc < 0.05) %>% 
        dplyr::arrange(p_val) %>% 
        group_by(mid, add = TRUE) %>%
        dplyr::slice(seq_len(n_dd[.$sid[1], .$i[1]])) %>% 
        summarize(hit = list(hit))) %>% 
    group_by(mid) %>% 
    summarize(hit = list(purrr::reduce(hit, c))) %>% 
    mutate_at("mid", as.character)

df <- fromList(set_names(top$hit, top$mid)) %>% 
    dplyr::mutate(
        code = apply(.[top$mid], 1, paste, collapse = ""),
        degree = apply(.[top$mid], 1, sum),
        hit = unique(unlist(top$hit))
    ) %>% {
        m <- match(.$hit, res$hit)
        dplyr::mutate(., 
            sid = res$sid[m], 
            i = res$i[m], 
            cat = res$category[m])
    } %>% 
    add_count(code) %>% 
    group_by(code) %>% 
    ungroup

m <- match(unique(df$code), df$code)
keep <- pull(top_n(df[m, ], 40, n), "code")
df <- dplyr::filter(df, code %in% keep)
m <- match(unique(df$code), df$code)
o <- order(df$degree[m], -df$n[m])

max <- ceiling(max(df$n)/1e3)*1e3
p1 <- ggplot(df, aes(x = code)) +
    stat_count(aes(fill = cat)) + 
    scale_x_discrete(limits = df$code[m][o], expand = c(0,0)) +
    scale_y_continuous(limits = c(0, max), expand = c(0,0)) +
    coord_trans(y = "sqrt", clip = "off") +
    scale_fill_manual(NULL, values = .cat_cols, 
        labels = c("EE", "DE", "DP", "DM", "DB")) +
    .prettify("classic") + theme(
        legend.margin = margin(0,0,0,0, "mm"),
        plot.margin = unit(c(2,2,2,2), "mm"),
        legend.key.size = unit(4, "mm"),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = NULL)

lgd <- get_legend(p1 + theme(legend.margin = margin(0,0,0,0, "mm")))
p1 <- p1 + theme(legend.position = "none")

dfm <- melt(df, variable.name = "method",
    id.var = setdiff(names(df), top$mid)) %>% 
    group_by(code, method) %>% 
    dplyr::slice(1)
p2 <- ggplot(dfm, aes(x = code, y = method, color = factor(value))) +
    scale_x_discrete(limits = df$code[m][o]) +
    scale_y_discrete(limits = rev(levels(res$mid))) +
    scale_color_manual(values = c("0" = "grey90", "1" = "black")) +
    geom_point(shape = 16, size = 1) +
    geom_path(size = 0.2, data = dplyr::filter(dfm, value != 0), aes(group = code)) +
    annotate("rect", alpha = 0.08, xmin = 0.5, xmax = Inf, 
        ymin = seq(0.5,nrow(top),2), ymax = seq(1.5, nrow(top)+1, 2)) +
    .prettify("classic") + theme(
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = NULL) 

anno <- read.csv(config$mids)
anno$type <- factor(anno$type, 
    levels = names(.typ_cols))

p3 <- ggplot(anno, aes(x = 0, y = id, fill = type)) +
    scale_fill_manual("method", values = .typ_cols, labels = .typ_labs) +
    geom_tile(col = "white") + coord_fixed(0.6) +
    scale_y_discrete(limits = rev(top$mid)) +
    .prettify("classic") + theme(aspect.ratio = NULL,
        legend.margin = margin(0,-4,0,0, "mm"),
        plot.margin = unit(c(0,6,0,0), "mm"),
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
        align = "hv", axis = "tl", 
        rel_heights = c(3, 2)),
    nrow = 1, rel_widths = c(7, 1))

ggsave(args$fig, p, dpi = 300,
    width = 15, height = 8, units = "cm")

