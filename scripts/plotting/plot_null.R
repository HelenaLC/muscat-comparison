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

# tidy results
#res <- list.files("results/kang/nill/", full.names = TRUE)
pvs <- .tidy_res(snakemake@input$res) %>% 
    setDT %>% split(
        by = c("sim_rep", "method"),  
        sorted = TRUE, flatten = FALSE) %>% 
    map_depth(2, "p_val") %>% 
    unlist

df <- data.frame(
    row.names = NULL,
    stringsAsFactors = FALSE,
    replicate = factor(gsub("(^[0-9]).*", "\\1", names(pvs))),
    method = gsub("[0-9]\\.([[:alpha:][:punct:]]+)[0-9]+", "\\1", names(pvs)),
    p_val = pvs) %>% 
    mutate_at("method", factor, levels = names(method_colors)) %>% 
    mutate_at("method", droplevels)

# plotting
p <- ggplot(df, aes(x = p_val, y = ..ndensity.., 
    col = method, fill = method, lty = replicate)) +
    geom_density(adjust = 0.2, alpha = 0.2, size = 0.3) + 
    facet_wrap(~ method, nrow = 3) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    guides(col = FALSE, lty = guide_legend(ncol = 1),
        fill = guide_legend(ncol = 3, order = 1,
            override.aes = list(alpha = 1, col = NA))) + 
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0.04)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0, 0.06)) +
    labs(x = "p-value", y = "normalized density") +
    prettify(theme = "classic",
        aspect.ratio = 2/3,
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom",
        legend.direction = "horizontal") +
    theme(panel.spacing = unit(0.2, "cm"))

# save data.frame to .rds & figure to .pdf
saveRDS(df, snakemake@output$res)
ggsave(snakemake@output$fig, p, 
    width = 15, height = 7.5, units = "cm",
    dpi = 300, useDingbats = FALSE)
