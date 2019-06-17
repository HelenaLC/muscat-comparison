# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(ggplot2)
        library(purrr)
    }))

# source utils
source(snakemake@input$utils)

# load truth & tidy results (split by gene category)
sim <- snakemake@input$sim#list.files("data/sim_data/kang", "d[a-z][0-9]+;", full.names = TRUE)
res <- snakemake@input$res#lapply(grep("d[a-z][0-9]+$", list.dirs("results/kang", full.names = TRUE), value = TRUE), list.files, full.names = TRUE) %>% unlist

cats <- c("de", "dp", "dm", "db")
names(cats) <- cats
gis <- lapply(cats, function(c)
    .get_gis(grep(paste0(c, "[0-9]+;"), sim, value = TRUE)))
dts <- lapply(cats, function(c)
    .tidy_res(grep(c, res, value = TRUE), gis[[c]]) %>% 
        setDT %>% split(
            sorted = TRUE, flatten = FALSE,
            by = c("sim_rep", "method")))

# caluclate performance for each
# p-adjustment type, method, and replicate 
ps <- paste0("p_adj.", c("loc", "glb"))
names(ps) <- ps
perf <- lapply(cats, function(c) {
    lapply(ps, function(p) {
        ps <- c("p_val", p)
        names(ps) <- c("p_val", "p_adj")
        lapply(seq_along(dts[[c]]), function(i)
            lapply(ps, function(p)
                dts[[c]][[i]] %>% 
                    map_depth(1, pull, p) %>% 
                    data.frame(check.names = FALSE)) %>% 
                c(list(truth = data.frame(
                    is_de = gis[[c]][[i]]$is_de)))) %>% 
            map(.calc_perf)
    })
})

# prep. data.frame for plotting
df <- map_depth(perf, 3, function(u)
    u %>% fdrtpr %>%
        select(method, thr, TPR, FDR) %>% 
        mutate_if(is.factor, as.character)) %>% 
    map_depth(2, bind_rows, .id = "rep") %>% 
    map_depth(1, bind_rows, .id = "p_adj") %>% 
    bind_rows(.id = "cat") %>% 
    # average over replicates
    group_by(cat, method, thr, p_adj) %>% 
    summarize_at(c("TPR", "FDR"), mean) %>% ungroup %>% 
    mutate_at("thr", function(u)
        as.numeric(gsub("thr", "", u))) %>%
    mutate_at("method", factor, levels = names(method_colors)) %>% 
    mutate_at("cat", factor, levels = cats, labels = toupper(cats)) %>% 
    mutate_at("p_adj", factor, levels = ps)

p <- .plot_perf_points(df, facet = c("p_adj", "cat"))
p$facet$params$ncol <- 4

saveRDS(p, snakemake@output$gg)
ggsave(snakemake@output$fig, p,
    width = 15, height = 12, units = "cm",
    dpi = 300, useDingbats = FALSE)
