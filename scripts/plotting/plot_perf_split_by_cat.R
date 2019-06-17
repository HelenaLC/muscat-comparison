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

# load truth & tidy results (split by gene category)
sim <- snakemake@input$sim#list.files("data/sim_data/kang", "d[a-z][0-9]+;", full.names = TRUE)
res <- snakemake@input$res#grep("d[a-z][0-9]+$", list.dirs("results/kang", full.names = TRUE), value = TRUE)

cats <- c("de", "dp", "dm", "db")
names(cats) <- cats
gis <- lapply(cats, function(c)
    .get_gis(grep(paste0(c, "[0-9]+;"), sim, value = TRUE)))
dts <- lapply(cats, function(c)
    .tidy_res(grep(c, res, value = TRUE), gis[[c]]) %>% 
        setDT %>% split(
            sorted = TRUE, flatten = FALSE,
            by = c("sim_rep", "method")))

# calculate performance across methods
# for ea. gene category & simulation replicate
ps <- c("p_val", "p_adj.loc")
names(ps) <- c("p_val", "p_adj")
res <- lapply(cats, function(c)
    lapply(seq_along(dts[[c]]), function(i)
        lapply(ps, function(p) 
            dts[[c]][[i]] %>% map(pull, p) %>% 
                data.frame(check.names = FALSE)) %>% 
            c(list(truth = data.frame(
                stringsAsFactors = FALSE,
                is_de = gis[[c]][[i]]$is_de,
                sim_rep = i, cat = c)))))

# prep. data.frame for plotting
df <- map_depth(res, 2, function(u)
    .calc_perf(u) %>% fdrtpr %>% 
        select(method, thr, TPR, FDR)) %>% 
    map(bind_rows, .id = "rep") %>% 
    bind_rows(.id = "cat") %>% 
    group_by(cat, method, thr) %>% 
    summarize_at(c("FDR", "TPR"), mean) %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>% 
    mutate_at("method", factor, levels = names(method_colors)) %>% 
    mutate_at("cat", factor, levels = cats, labels = toupper(cats))

p <- .plot_perf_points(df, facet = "cat")
p$facet$params$ncol <- nlevels(df$cat)
     
saveRDS(p, snakemake@output$gg)
ggsave(snakemake@output$fig, p,
    width = 15, height = 8, units = "cm",
    dpi = 300, useDingbats = FALSE)
