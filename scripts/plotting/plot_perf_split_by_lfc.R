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
gis <- .get_gis(snakemake@input$sim)#list.files("data/sim_data/kang", "de10;", full.names = TRUE))
dt <- .tidy_res(snakemake@input$res, gis) %>%#.tidy_res(list.files("results/kang/de10", full.names = TRUE), gis) %>% 
    setDT %>% split(
        sorted = TRUE, flatten = FALSE,
        by = c("sim_rep", "method"))

# specify logFC-groups
fil <- list("abs(u)<=1", "abs(u)>1&abs(u)<=2", "abs(u)>2")
names(fil) <- c("low", "mid", "high")

# calculate performance
ps <- c("p_val", "p_adj.loc")
names(ps) <- c("p_val", "p_adj")
res <- lapply(seq_along(dt), function(i) {
    re <- lapply(seq_along(fil), function(j) {
        re <- lapply(ps, function(p) 
            dt[[i]] %>% map(pull, p) %>% bind_cols %>% 
                data.frame(check.names = FALSE))
        u <- gis[[i]]$sim_logFC
        keep <- eval(parse(text = fil[j]))
        keep[is.na(keep)] <- TRUE
        lapply(re, function(u) u[keep, ]) %>% 
            c(list(truth = data.frame(
                is_de = gis[[i]]$is_de[keep],
                group = names(fil)[j],
                stringsAsFactors = FALSE)))
    })
    lapply(names(re[[1]]), function(i)
        map(re, i) %>% bind_rows) %>% 
        set_names(c(ps, "truth"))
})

# calculate performance
perf <- map(res, function(u)
    .calc_perf(u, facet = "group") %>% fdrtpr %>% 
        mutate(group = gsub("group:", "", splitval)) %>% 
        select(c("method", "group", "thr", "FDR", "TPR")))

# prep. data.frame for plotting
df <- perf %>% 
    bind_rows(.id = "rep") %>% 
    group_by(group, method, thr) %>% 
    summarize_at(c("FDR", "TPR"), mean) %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    mutate_at("thr", function(u) as.numeric(gsub("thr", "", u))) %>% 
    mutate_at("method", factor, levels = names(method_colors)) %>% 
    mutate_at("group", factor, 
        levels = c("low", "mid", "high", "overall"),
        labels = c("|logFC| <= 1", "1 < |logFC| <= 2", "|logFC| > 2", "overall"))

p <- .plot_perf_points(df, facet = "group")
p$facet$params$ncol <- nlevels(df$group)

ggsave(snakemake@output$fig, p,
    width = 15, height = 7, units = "cm",
    dpi = 300, useDingbats = FALSE)
