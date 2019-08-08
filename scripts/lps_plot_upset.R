source(snakemake@config$utils)

# load packages
suppressMessages({
    library(dplyr)
    library(purrr)
    library(UpSetR)
})

# load results
snakemake@input$res

fns <- list.files("results/lps", full.names = TRUE)
res <- lapply(fns, readRDS) %>% map("tbl") %>% 
    set_names(gsub("(.rds)", "", basename(fns)))

rmv <- vapply(res, inherits, what = "error", logical(1))
df <- bind_rows(res[!rmv]) %>% 
    mutate_at("method_id", factor, levels = names(.meth_cols)) %>% 
    mutate_at("method_id", droplevels)

gs <- filter(df, p_adj.loc < 0.05) %>% 
    split(.$method_id) %>% 
    map(pull, "gene")

p <- upset(fromList(gs), 
    sets = levels(df$method_id),
    nintersects = 80,
    mb.ratio = c(0.6, 0.4))
p    
    
