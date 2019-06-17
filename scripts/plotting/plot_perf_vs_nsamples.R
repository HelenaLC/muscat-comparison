# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(iCOBRA)
        library(ggplot2)
        library(magrittr)
        library(purrr)
        library(SingleCellExperiment)
    })
)

# source utils
source(snakemake@input$utils)

# load truth & tidy results
#sim <- list.files("data/sim_data/kang", "ns", full.names = TRUE)
#res <- list.files("results/kang/de10_ns", full.names = TRUE)
gis <- .get_gis(snakemake@input$sim)#.get_gis(sim)#
dt <- .tidy_res(snakemake@input$res, gis) %>% #.tidy_res(res, gis) %>% #
    setDT %>% split(
        sorted = TRUE, flatten = FALSE,
        by = c("sim_rep", "method", "n_samples"))

# calculate performance
ps <- c("p_val", "p_adj.loc")
names(ps) <- c("p_val", "p_adj")
methods <- names(dt[[1]])
perf <- lapply(seq_along(dt), function(i) {
    is_de <- gis[[i]]$is_de
    lapply(ps, function(p)
        dt[[i]] %>% map_depth(2, pull, p) %>% 
            map_depth(1, data.frame, check.names = FALSE) %>% 
            bind_rows) %>% 
        c(list(truth = data.frame(
            stringsAsFactors = FALSE,
            is_de = rep(is_de, length(methods)), 
            method = rep(methods, each = length(is_de)))))
}) %>% 
    map_depth(2, data.frame, check.names = FALSE) %>% 
    map(.calc_perf, facet = "method")

# prep. data.frame for plotting
df <- map(perf, fdrtpr) %>% 
    map(select, method, splitval, FDR, TPR, thr) %>% 
    bind_rows(.id = "rep") %>% 
    mutate(n_samples = method, method = splitval) %>% 
    filter(method != "overall") %>% 
    mutate_at("method", function(u) 
        factor(gsub("^method:", "", u), 
            levels = names(method_colors))) %>% 
    mutate_at("thr", function(u) 
        as.numeric(gsub("thr", "", u))) %>%
    mutate_at("n_samples", function(u) factor(u, 
        levels = sort(as.numeric(unique(u))))) %>% 
    group_by(thr, method, n_samples) %>% 
    summarise_at(c("FDR", "TPR"), mean)

# ------------------------------------------------------------------------------
# TPR vs. FDR points at thresholds 0.01, 0.05, and 0.1 
# ------------------------------------------------------------------------------
#   - facetted by nb. of samples
#   - colored by method
#   - averaged across replicate
# ------------------------------------------------------------------------------

p <- .plot_perf_points(df, facet = "n_samples")
p$facet$params$ncol <- nlevels(df$n_samples)
saveRDS(p, snakemake@output$gg)
ggsave(snakemake@output$fig1, p,
    width = 15, height = 8, units = "cm", 
    dpi = 300, useDingbats = FALSE) 

# ------------------------------------------------------------------------------
# TPR vs. FDR curves 
# ------------------------------------------------------------------------------
#   - facetted by method
#   - colored by nb. of samples
#   - one curve per replicate
# ------------------------------------------------------------------------------

lvls <- sort(as.numeric(as.character(levels(df$n_samples))))
cols <- setNames(brewer.pal(nlevels(df$n_samples), "RdYlBu"), lvls)

u <- perf[[1]]
u@fdrtpr <- bind_rows(map(perf, function(u) u@fdrtpr), .id = "rep")
u@fdrtprcurve <- bind_rows(map(perf, function(u) u@fdrtprcurve), .id = "rep")

u@fdrtpr <- mutate_at(u@fdrtpr,
    c("splitval", "fullmethod"), 
    function(u) factor(
        gsub("method:", "", u), 
        levels = names(method_colors)))
u@fdrtprcurve <- mutate_at(u@fdrtprcurve, 
    c("splitval", "fullmethod"), 
    function(u) factor(
        gsub("method:", "", u),
        levels = names(method_colors)))

u <- u %>% reorder_levels(lvls)
u <- prepare_data_for_plot(u, colorscheme = cols, 
    incloverall = FALSE, facetted = TRUE)
p <- .plot_perf_curves(u) + theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.text = element_text(size = 4),
    panel.spacing = unit(0.2, "cm"))
p$facet$params$ncol <- 6
p$layers[[1]]$aes_params$size   <- 0.2 # dashed lines at thrs
p$layers[[2]]$aes_params$size   <- 0.2 # line size
p$layers[[2]]$aes_params$alpha  <- 0.8 # line alpha
p$layers[[2]]$show.legend <- FALSE
p$layers[[3]]$aes_params$stroke <- 0.1 # point size
p$layers[[3]]$aes_params$alpha  <- 0.8 # point alpha
p$guides$colour <- guide_legend("n_samples", override.aes = list(alpha = 1, size = 2))

ggsave(snakemake@output$fig2, p,
    width = 15, height = 14, units = "cm", 
    dpi = 300, useDingbats = FALSE)
