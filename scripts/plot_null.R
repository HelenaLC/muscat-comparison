suppressMessages({
    library(data.table)
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
    library(reshape2)
})

# id <- "magl"
# args <- list(
#     res = list.files("results", sprintf("%s,nill", id), full.names = TRUE),
#     ggp = file.path("plots", sprintf("%s-null.rds", id)),
#     fig = file.path("plots", sprintf("%s-null.pdf", id)))

df <- .read_res(args$res) %>% 
    dplyr::filter(!is.na(p_val)) %>% 
    dplyr::rename(method = mid, replicate = i)

p <- ggplot(df, aes(x = p_val, y = ..ndensity.., 
    col = method, fill = method, lty = replicate)) +
    facet_wrap(~ method, ncol = 4) + 
    geom_density(adjust = 0.2, size = 0.3, alpha = 0.1) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    scale_color_manual(values = .meth_cols) +
    scale_fill_manual(values = .meth_cols) +
    guides(col = FALSE,
        lty = guide_legend(ncol = 1, order = 2),
        fill = guide_legend(ncol = 1, order = 1,
            override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("normalized density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + theme(aspect.ratio = 1/2,
        #legend.position = "bottom",
        legend.box.just = "left",
        panel.grid = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(color = "grey"),
        strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

saveRDS(p, args$ggp)
ggsave(args$fig, p,
    units = "cm", width = 15, height = 8,
    dpi = 300, useDingbats = FALSE)

