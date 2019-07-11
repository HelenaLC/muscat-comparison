source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
    library(reshape2)
})

#fns <- list.files("~/documents/kang", "nill", full.names = TRUE)
df <- .read_res(snakemake@input$res) %>% 
    dplyr::filter(!is.na(p_val)) %>% 
    dplyr::rename(method = mid, replicate = i)

p <- ggplot(df, aes(x = p_val, y = ..ndensity.., 
    col = method, fill = method, lty = replicate)) +
    facet_wrap(~ method, nrow = 3) + 
    geom_density(adjust = 0.2, size = 0.3, alpha = 0.1) +
    scale_color_manual(values = .meth_cols) +
    scale_fill_manual(values = .meth_cols) +
    guides(col = FALSE,
        lty = guide_legend(ncol = 1, order = 2),
        fill = guide_legend(ncol = 3, order = 1,
            override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + theme(aspect.ratio = 1/2,
        legend.position = "bottom",
        legend.box.just = "left",
        panel.grid = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(color = "grey"),
        strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

saveRDS(p, snakemake@output$ggp)
ggsave(snakemake@output$fig, p,
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)

