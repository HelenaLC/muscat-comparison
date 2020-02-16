suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(purrr)
})

lapply(config$dids, function(id) {
    fns <- sprintf("%s,all-perf_by_cat_%s.rds", id, c("loc", "glb"))
    ps <- lapply(file.path("plots", fns), readRDS)
    df <- bind_rows(map(ps, "data"), .id = "p_adj") %>% 
        mutate_at("p_adj", factor, labels = paste0("p_adj.", c("loc", "glb")))
    p <- .plot_perf_points(df, facet = c("p_adj", "splitval"))
    p$facet$params$ncol <- 4
    fns <- sprintf("%s-perf_by_cat%s", id, c(".rds", ".pdf"))
    fns <- file.path("figures", fns)
    saveRDS(p, fns[1])
    ggsave(fns[2], p, 
        dpi = 300, useDingbats = FALSE,
        width = 15, height = 10, units = "cm")
})