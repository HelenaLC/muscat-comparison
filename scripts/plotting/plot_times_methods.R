# source utils
source(snakemake@input$dir_utils)

library(ggplot2)
df <- read.csv(snakemake@input$dir_res)
df$id <- with(df, paste0(n_genes, method))
p <- ggplot(df, aes(x = n_genes, y = time, group = id, col = method)) +
    geom_boxplot() + 
    scale_x_continuous(expand = c(.02,0), trans = "log2", 
        breaks = unique(df$n_genes), labels = unique(df$n_genes)) +
    scale_y_continuous(expand = c(.03,0), trans = "log2", 
        breaks = 2^log2(c(1, 10, 60, 180))) +
    labs(y = "time (s)") + prettify + theme(aspect.ratio = 2/3) + 
    stat_smooth(formula = y ~ x, aes(group = method), alpha = 0.1, size = 0.2)
ggsave(snakemake@output$dir_fig, p, width = 18, height = 10, units = "cm")
