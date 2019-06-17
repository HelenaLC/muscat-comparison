# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(ggplot2)
    }))

# source utils
source(snakemake@input$dir_utils)

# load results
times <- read.csv(snakemake@input$dir_res)

df <- times
df$id <- with(df, paste0(n_cells, method))
p <- ggplot(df, aes(x = n_cells, y = time, group = id, col = method)) +
    geom_boxplot() + 
    scale_x_continuous(expand = c(.02,0), trans = "log2", 
        breaks = unique(df$n_cells), labels = unique(df$n_cells) * 12) +
    scale_y_continuous("time (s)", expand = c(.03,0)) +
    prettify() + theme(aspect.ratio = 2/3) + 
    stat_smooth(formula = y ~ x, aes(group = method), alpha = 0.1, size = 0.2)
ggsave(snakemake@output$dir_fig, p, width = 18, height = 10, units = "cm")
