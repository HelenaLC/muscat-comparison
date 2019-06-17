source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(dplyr)
    library(iCOBRA)
    library(ggplot2)
    library(purrr)
})

res <- lapply(snakemake@input$res, readRDS) %>%
    map("tbl") %>% bind_rows %>% split(.$mid)

cd <- COBRAData( 
    pval = as.data.frame(bind_cols(map(res, "p_val"))),
    padj = as.data.frame(bind_cols(map(res, "p_adj.loc"))),
    truth = data.frame(is_de = res[[1]]$is_de))

perf <- calculate_performance(
    cobradata = cd, 
    binary_truth = "is_de", 
    aspects = c("fdrtpr", "fdrtprcurve"))

p <- plot_fdrtprcurve(
    prepare_data_for_plot(perf), 
    linewidth = 0.6, pointsize = 1) +
    guides(color = guide_legend(override.aes = list(size = 1.2, linetype = 0))) +
    scale_color_manual(values = .meth_cols) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0.05, 0)) +
    scale_x_sqrt(breaks = c(0.01, 0.05, 0.1, seq(0.2, 1, 0.2)), expand = c(0.05, 0)) + 
    .prettify("bw", panel.grid = element_blank()) + theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
p$layers[[1]] <- NULL # remove vertical dashed lightgrey lines 
p$layers[[1]]$aes_params$size <- 0.4 # dashed lines at thresholds

ggsave(snakemake@output$fig, p,
    units = "cm", width = 12, height = 8,
    dpi = 200, useDingbats = FALSE)
