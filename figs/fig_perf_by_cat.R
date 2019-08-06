config <- yaml::read_yaml("config.yaml")
source(config$utils)

suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(purrr)
})

fns <- sprintf("perf_by_cat_%s.rds", c("loc", "glb"))
for (id in config$dids) {
    ps <- lapply(file.path("figures", id, fns), readRDS)
    df <- bind_rows(map(ps, "data"), .id = "p_adj") %>% 
        mutate_at("p_adj", factor, labels = paste0("p_adj.", c("loc", "glb")))
    p <- .plot_perf_points(df, facet = c("p_adj", "splitval"))
    p$facet$params$ncol <- 4
    
    saveRDS(p, file.path("figures", id, "perf_by_cat.rds"))
    ggsave(file.path("figures", id, "perf_by_cat.pdf"), 
        p, dpi = 300, useDingbats = FALSE,
        width = 15, height = 10.1, units = "cm")
}






