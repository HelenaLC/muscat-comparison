# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(ggplot2)
        library(iCOBRA)
        library(purrr)
    })
)

# source utils
source(snakemake@input$utils)

# load truth & tidy results
gis <- .get_gis(snakemake@input$sim)
dt <- .tidy_res(snakemake@input$res, gis) %>% 
    setDT %>% split(
        sorted = TRUE, flatten = FALSE,
        by = "method")

# prep. data for iCOBRA
is_de <- gis[[1]]$is_de
ps <- c("p_val", "p_adj.loc")
names(ps) <- c("p_val", "p_adj")
res <- lapply(ps, function(p)
    dt %>% map(pull, p) %>% bind_cols %>% 
    data.frame(check.names = FALSE)) %>% 
    c(list(truth = data.frame(is_de)))

# calculate performance & plot FDR vs. TPR curve
p <- res %>% .calc_perf %>% 
    prepare_data_for_plot(colorscheme = method_colors) %>% 
    reorder_levels(names(method_colors)) %>% 
    .plot_perf_curves + theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +
    guides(color = guide_legend(ncol = 3))

ggsave(snakemake@output$fig, p,
    width = 15, height = 8, units = "cm",
    dpi = 300, useDingbats = FALSE)