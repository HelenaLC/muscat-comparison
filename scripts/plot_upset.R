source(snakemake@config$utils)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(purrr)
})

fns <- list.files("results/kang", "d[a-z]10;", full.names = TRUE)
res <- lapply(fns, readRDS) %>% map("tbl") %>% bind_rows %>% 
    dplyr::mutate(hit = paste(gene, cluster_id, sid, i, sep = ";"))

n_de <- res %>% 
    dplyr::filter(mid == .$mid[1]) %>% 
    group_by(sid, i) %>% 
    summarize(n_de = sum(is_de)) %>% 
    acast(sid ~ i, value.var = "n_de")

top <- res %>% group_by(sid, i) %>% do({
    n <- n_de[.$sid[1], .$i[1]]
    group_by(., mid, add = TRUE) %>%
        arrange(p_val) %>%
        slice(seq_len(n)) %>% 
        summarize(hit = list(hit))
}) %>% group_by(mid) %>% summarize(hit = list(reduce(hit, c)))

df <- fromList(set_names(top$hit, top$mid)) %>% 
    mutate(hit = unique(unlist(top$hit))) %>% 
    mutate(cat = res$category[match(hit, res$hit)])

fun <- function(row, cat) {
    data <- row["cat"] == cat
}
cols <- c("royalblue", "cornflowerblue", "red3", "tomato", "orange", "gold")
names(cols) <- levels(res$category)
upset(df, sets = top$mid[1:5], nintersects = 20, query.legend = "top",
    queries = lapply(levels(droplevels(res$category)), function(c)
        list(query = fun, params = list(c), color = cols[c], active = TRUE, query.name = c)))



 