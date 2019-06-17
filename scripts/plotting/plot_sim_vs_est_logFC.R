# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(cowplot)
        library(dplyr)
        library(ggplot2)
        library(purrr)
    }))

# source utils
source(snakemake@input$utils)

# load truth & tidy results (split by gene category)
sim <- snakemake@input$sim#list.files("data/sim_data/kang", "d[a-z][0-9]+;", full.names = TRUE)
#res <- grep("d[a-z][0-9]+$", list.dirs("results/kang", full.names = TRUE), value = TRUE)
res <- snakemake@input$res#lapply(res, list.files, full.names = TRUE) %>% unlist

cats <- c("de", "dp", "dm", "db")
names(cats) <- cats
gis <- lapply(cats, function(c)
    .get_gis(grep(paste0(c, "[0-9]+;"), sim, value = TRUE)))
dts <- lapply(cats, function(c)
    .tidy_res(grep(c, res, value = TRUE), gis[[c]])) %>% 
    bind_rows(.id = "cat")

# remove methods w/o any logFC estimates
rmv <- dts %>% 
    group_by(method) %>% 
    group_split() %>% 
    set_names(levels(dts$method)) %>% 
    map(function(u) all(is.na(u$est_logFC))) %>% 
    unlist %>% (function(u) names(u)[u])

res <- dts %>%
    filter(!method %in% rmv) %>% 
    mutate_at("method", droplevels) 
    
# compute residual sum-of-squares
# b/w simulated & estimated logFCs
r2 <- res %>% 
    filter(is_de == 1) %>% 
    group_by(cat, method) %>% 
    summarize(r2 = sum((est_logFC - sim_logFC) ^ 2)) %>% 
    split(.$cat)

# prep. data.frame for plotting
df <- res %>% 
    mutate_at("is_de", as.factor) %>% 
    mutate_at("est_logFC", function(u) { 
        u[u > 6] <- 6; u[u < -6] <- -6; u })

# sample 10k points per method
df <- df %>% group_by(cat, method) %>% 
    group_split() %>% map(function(u) 
        u[sample(nrow(u), 5e3), ]) %>% 
    bind_rows %>% split(.$cat)

# plot scatter of sim_logFC vs. est_logFC,
# facetted by method, colored by is_de
cats <- c("de", "dp", "dm", "db")
names(cats) <- cats
ps <- lapply(cats, function(c) {
    ggplot(df[[c]], aes(x = est_logFC, y = sim_logFC, 
        col = is_de, size = is_de, alpha = is_de)) +
        geom_point() +
        geom_abline(intercept = 0, slope = 1, size = 0.1) +
        geom_text(data = r2[[c]], size = 1.5, inherit.aes = FALSE,
            x = -6.5, y = 6.5, hjust = 0, vjust = 1, 
            aes(label = round(r2, 2))) +
        facet_wrap(~method, ncol = 7) +
        scale_size_manual(values = c("0" = 0.01, "1" = 0.2)) +
        scale_alpha_manual(values = c("0" = 0.01, "1" = 0.2)) +
        scale_color_manual(values = c("0" = "royalblue", "1" = "tomato")) +
        guides(color = guide_legend("is_de",
            override.aes = list(alpha = 1, size = 2))) +
        scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, 3), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-7, 7), breaks = seq(-6, 6, 3), expand = c(0, 0)) +
        prettify(theme = "bw", aspect.ratio = 1,
            legend.position = "bottom", 
            legend.direction = "horizontal") +
        theme(strip.text = element_text(size = 3))
})

# arrange panels
lgd <- get_legend(ps[[1]])
ps <- lapply(ps, function(p) p + theme(legend.position = "none"))
ps <- c(lapply(ps[-4], function(p) p + labs(x = NULL)), ps[4])
p <- plot_grid(axis = "t", align = "v",
    plotlist = c(ps, list(lgd)), ncol = 1, 
    labels = toupper(cats), label_size = 8,
    rel_heights = c(5,5,5,5.3,0.3))

ggsave(snakemake@output$fig, p,
    width = 12, height = 20, units = "cm",
    dpi = 300, useDingbats = FALSE)
