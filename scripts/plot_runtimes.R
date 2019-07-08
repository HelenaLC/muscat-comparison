source(snakemake@config$utils)

suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(purrr)
})

fns <- list.files("/Users/helena/Dropbox/portmac/results/kang", "ds10_nc;", full.names = TRUE)
res <- lapply(fns, readRDS)
tbl <- map(res, "tbl")
rmv <- vapply(tbl, inherits, what = "error", logical(1))
nms <- basename(fns)[!rmv]

rts <- map(res[!rmv], "rt") %>% 
    vapply(sum, numeric(1)) %>% 
    set_names(nms)

df <- map(tbl[!rmv], mutate_if, is.factor, as.character) %>% 
    set_names(nms) %>% bind_rows(.id = "id")

df <- df[match(nms, df$id), ] %>% 
    mutate(t = rts) %>% 
    mutate(id = paste(mid, c, sep = "--")) %>% 
    mutate_at("c", as.numeric)

ggplot(df, aes(x = c, y = t, group = id, col = mid, fill = mid)) +
    geom_smooth(method = "lm", aes(group = mid), size = 0.2, alpha = 0.05) + 
    geom_boxplot(alpha = 0.8) + scale_y_log10() +
    scale_x_log10(breaks = unique(df$c)) +
    scale_color_manual(NULL, values = .meth_cols) +
    scale_fill_manual(NULL, values = .meth_cols, guide = FALSE) +
    labs(x = "#(cells) per cluster-sample", y = "runtime (s)") +
    .prettify("bw") + theme(aspect.ratio = 2/3)

ggsave(snakemake@output$fig, p,
    units = "cm", width = 7, height = 4,
    dpi = 300, useDingbats = FALSE)    
