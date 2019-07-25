source(snakemake@config$utils)

suppressMessages({
    library(cowplot)
    library(dplyr)
    library(ggplot2)
    library(purrr)
})

#fns <- list.files("results/kang", "ds10_n[g|c];", full.names = TRUE)
fns <- snakemake@input$res
pat <- "%s;%s;%s;%s;g%s;c%s;k%s;s%s"
tbl <- .read_res(fns) %>% 
    dplyr::mutate(id = sprintf(pat, sid, i, mid, j, g, c, k, s))
rts <- .read_res(fns, slot = "rt") %>% 
    vapply(sum, numeric(1)) %>% 
    set_names(gsub(".rds", "", basename(fns)))

m <- match(tbl$id, names(rts))
df <- data.frame(tbl, rt = rts[m],
    stringsAsFactors = FALSE) %>% 
    mutate_at(c("c", "g"), function(u)
        as.numeric(as.character(u))) %>% 
    mutate_at("c", function(u) u * 2 * 2 * 3) %>% 
    mutate_at("g", function(u) u * 2)

ps <- lapply(c("c", "g"), function(x) {
   u <- dplyr::filter(df, sid == paste0("ds10_n", x)) %>% 
       dplyr::mutate(id = paste(mid, get(x), sep = "--")) %>% 
       group_by(id, i) %>% dplyr::slice(1)
    
    ss <- strsplit(unique(u$id), "--")
    m <- sapply(ss, .subset, 1)
    n <- sapply(ss, .subset, 2)
    o1 <- order(as.numeric(n))
    o2 <- match(m, levels(u$mid))
    lvls <- unique(u$id)[order(o2)[o1]]
    u$id <- factor(u$id, levels = lvls)
    
    ggplot(u, aes_string(x = x, y = "rt", group = "id", col = "mid", fill = "mid")) +
        geom_boxplot(alpha = 0.8, size = 0.4, outlier.size = 0.4) + 
        scale_y_log10() + 
        scale_x_log10(breaks = unique(u[[x]])) +
        scale_color_manual(NULL, values = .meth_cols) +
        scale_fill_manual(NULL, values = .meth_cols) +
        guides(col = guide_legend(ncol = 3)) +
        labs(x = sprintf("nb. of %s", c(c = "cells", g = "genes")[x]), y = "runtime (s)") +
        .prettify("bw") + theme(legend.position = "bottom", aspect.ratio = 2/3)
})

lgd <- get_legend(ps[[1]])
ps <- lapply(ps, "+", theme(legend.position = "none"))

p <- plot_grid(
    plot_grid(nrow = 1,
        plotlist = ps,  
        labels = "auto", 
        label_size = 10, 
        label_fontface = "bold"),
    lgd, ncol = 1, 
    rel_heights = c(4, 1))

ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 6.4,
    dpi = 300, useDingbats = FALSE)    

