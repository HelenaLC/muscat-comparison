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

# load truth & tidy results
gis <- .get_gis(snakemake@input$sim)
df <- .tidy_res(snakemake@input$res, gis) %>% 
    setDT %>% split(
        sorted = TRUE, flatten = FALSE,
        by = c("sim_rep", "method")) 

# get top ranked genes
orders <- map_depth(df, 2, select, "p_adj.loc") %>% 
    map(function(u) bind_cols(u) %>% set_colnames(names(u))) %>% 
    map(apply, 2, order)

top <- do.call("rbind", lapply(seq_along(orders), function(i) {
    gs <- with(gis[[i]], paste(i, cluster_id, gene))
    apply(orders[[i]], 2, function(j) gs[j][seq_len(100)])
}))

# prep. for plottting
tib <- tibble(
    method = colnames(top),
    hit = split(top, col(top))) %>% 
    unnest %>% group_by(hit) %>%
    summarize(method = list(method))

p <- ggplot(tib, aes(x = method, y = ..count../sum(..count..))) +
    geom_bar(fill = "royalblue") +
    scale_x_upset(
        order_by = "degree",
        n_intersections = 80) +
    scale_y_continuous(labels = function(u) 
        scales::percent(u, accuracy = 1),
        limits = c(0, 0.2), breaks = seq(0, 1, 0.05)) +
    labs(x = NULL, y = NULL) + 
    prettify() + theme(panel.grid.major.x = element_blank()) +
    theme_combmatrix(
        combmatrix.panel.line.size = 0.1,
        combmatrix.panel.point.size = 0.2,
        combmatrix.label.text = element_text(size = 5, 
            color = method_colors[rev(colnames(top))]),
        combmatrix.label.height = unit(2.5, "cm"))
p$coordinates$levels <- colnames(top)

ggsave(snakemake@output$fig, p,
    width = 15, height = 7.5, units = "cm", 
    dpi = 300, useDingbats = FALSE)  
