# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(ggplot2)
        library(ggupset)
        library(magrittr)
        library(purrr)
        library(SingleCellExperiment)
        library(tidyr)
    })
)

# source utils
source(snakemake@input$utils)

# load truth & tidy results (split by gene category)
sim <- list.files("data/sim_data/kang", "d[a-z][0-9]+;", full.names = TRUE)
res <- grep("d[a-z][0-9]+$", list.dirs("results/kang", full.names = TRUE), value = TRUE)
res <- lapply(res, list.files, full.names = TRUE) %>% unlist

#sim <- snakemake@input$sim
#res <- snakemake@input$res

cats <- c("de", "dp", "dm", "db")
names(cats) <- cats
gis <- lapply(cats, function(c)
    .get_gis(grep(paste0(c, "[0-9]+;"), sim, value = TRUE)))
dts <- lapply(cats, function(c)
    .tidy_res(grep(c, res, value = TRUE), gis[[c]]) %>% 
        setDT %>% split(
            sorted = TRUE, flatten = FALSE,
            by = c("sim_rep", "method")))

# get top ranked genes for ea. 
# gene category, method, and replicate
orders <- map_depth(dts, 3, select, "p_val") %>%
    map_depth(2, function(u) bind_cols(u) %>% set_colnames(names(u))) %>% 
    map_depth(2, apply, 2, order)

ns <- map_depth(gis, 2, function(u) sum(u$is_de))
top <- lapply(cats, function(c) {
    gs <- lapply(seq_along(orders[[c]]), function(i) {
        gs <- with(gis[[c]][[i]], paste(c, i, cluster_id, gene))
        apply(orders[[c]][[i]], 2, function(j) gs[j][seq_len(ns[[c]][[i]])])
    })
    do.call("rbind", gs)
})

# prep. tibble for plottting
is_de <- lapply(cats, function(c)
    lapply(names(gis[[c]]), function(i) {
        u <- gis[[c]][[i]]
        u$id <- paste(c, i, u$cluster_id, u$gene)
        return(u)
    }) %>% bind_rows) %>% bind_rows

tib <- lapply(top, function(u)
    tibble(
        method = colnames(u),
        hit = split(u, col(u))) %>% 
        unnest %>% group_by(hit) %>%
        summarize(method = list(method))) %>% 
    bind_rows(.id = "cat") %>% 
    mutate(is_de = is_de$is_de[match(hit, is_de$id)]) %>% 
    mutate_at("cat", factor, levels = cats, labels = toupper(cats))

p <- ggplot(tib, aes(x = method, y = ..count../sum(..count..))) +
    facet_wrap("cat", ncol = 1) +
    geom_bar(fill = "royalblue") +
    scale_x_upset(
        order_by = "degree",
        n_intersections = 50) +
    scale_y_continuous(labels = function(u) 
        scales::percent(u, accuracy = 1),
        limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
    labs(x = NULL, y = NULL) + 
    prettify(theme = "bw") + 
    theme(panel.grid.major.x = element_blank()) +
    theme_combmatrix(
        combmatrix.panel.line.size = 0.1,
        combmatrix.panel.point.size = 0.4,
        combmatrix.label.text = element_text(size = 6, 
            color = method_colors[rev(colnames(top[[1]]))]),
        combmatrix.label.height = unit(4, "cm"))
p$coordinates$levels <- colnames(top[[1]])

ggsave(snakemake@output$fig, p,
    width = 15, height = 16, units = "cm", 
    dpi = 300, useDingbats = FALSE)  

# sub <- filter(tib, cat == "DE")
# 
# # count occurences of ea. set
# mids <- names(method_colors)
# sets <- factor(sapply(sub$method, function(u) 
#     paste(u[match(mids, u, nomatch = 0)], collapse = ";")))
# n_de <- setNames(sapply(split(sub$is_de, sets), sum), levels(sets))
# set_counts <- table(sets)
# range(set_counts)
# 
# 
# # prep. data.frame for plotting
# sets_split <- strsplit(names(set_counts), split = ";")
# points <- t(sapply(sets_split, function(u) mids %in% u)) %>% 
#     set_colnames(mids) %>% `[`(TRUE, colSums(.) != 0)
# ms_keep <- colnames(points)
# df <- data.frame(
#     set = names(set_counts),
#     count = as.numeric(set_counts),
#     degree = sapply(sets_split, length),
#     points,
#     n_de = n_de[names(set_counts)],
#     check.names = FALSE)
# 
# df <- df[match(rev(levels(sets)), df$set, nomatch = 0), ]
# df <- arrange(df, desc(count))[seq_len(40), ]
# o <- df$set[order(df$degree)]
# df$set <- factor(as.character(df$set), levels = o)
# dfm <- melt(df, id.var = setdiff(colnames(df), mids))
# 
# foo <-  data.frame(x = 1:2, is_de = c("TRUE", "FALSE")) 
# foo <- ggplot(foo, aes(x, fill = is_de)) + 
#     stat_count() + prettify() +
#     scale_fill_manual(values = c("grey", "tomato"))
# lgd <- get_legend(foo)
# 
# p1 <- ggplot(df, aes(x = set)) +
#     geom_bar(stat = "identity", fill = "grey", aes(y = count/nrow(sub))) +
#     geom_bar(stat = "identity", fill = "tomato", aes(y = n_de/nrow(sub))) +
#     scale_x_discrete(limits = o) +
#     scale_y_continuous(labels = function(u)
#        scales::percent(u, accuracy = 1),
#         limits = c(0, 0.2), breaks = seq(0, 1, 0.05)) +
#     prettify(theme = "minimal") + theme(
#         axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         panel.grid.major.x = element_blank())
# 
# p2 <- ggplot(dfm, aes(x = set, y = variable, alpha = value)) +
#     scale_alpha_manual(values = c("FALSE" = 0.08, "TRUE" = 1)) +
#     annotate("rect", xmin = 0.5, xmax = Inf, 
#         ymin = seq(1.5, length(ms_keep), 2), 
#         ymax = seq(0.5, length(ms_keep), 2),
#         fill = "grey95") +
#     geom_point(size = 0.8) + geom_path(size = 0.1,
#         data = filter(dfm, value == 1), aes(group = set)) +
#     scale_x_discrete(limits = o) +
#     scale_y_discrete(limits = rev(ms_keep)) +
#     guides(alpha = FALSE) +
#     prettify(theme = "minimal") + theme(
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 5,
#             color = rev(method_colors[ms_keep])),
#         axis.title = element_blank(),
#         panel.grid = element_blank())
# 
# p <- plot_grid(
#     p1 + theme(plot.margin = unit(rep(0,4), "mm")), 
#     p2 + theme(plot.margin = unit(rep(0,4), "mm")), 
#     ncol = 1, axis = "l", align = "v")
# 
# ggsave("/users/helena/desktop/upset.pdf", p,
#     width = 15, height = 8, units = "cm", 
#     dpi = 300, useDingbats = FALSE) 

