options(conflicts.policy = list(warn = FALSE))

source(snakemake@config$utils)

suppressMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

wcs <- snakemake@wildcards
wcs <- wcs[names(wcs) != ""]
for (wc in names(wcs))
    assign(wc, wcs[[wc]])

sim <- readRDS(snakemake@input$sim)

meth_pars <- as.list(jsonlite::fromJSON(snakemake@input$meth_pars))
run_pars <- as.list(jsonlite::fromJSON(snakemake@input$run_pars))

set.seed(run_pars$seed + as.numeric(j))

# subset clusters & samples
kids <- levels(sim$cluster_id)
sids <- levels(sim$sample_id)
m <- match(sids, sim$sample_id)
gids <- sim$group_id[m]

if (k != "x") kids <- sample(kids, k)
if (s != "x") sids <- sapply(split(sids, gids), sample, s)

sim <- .filter_sce(sim, kids, sids)

# subset genes & cells
gs <- rownames(sim)
cs <- colnames(sim)

if (g != "x") 
    gs <- sample(gs, g)

if (c != "x") {
    cs <- split(cs, list(sim$cluster_id, sim$sample_id))
    cs <- unlist(lapply(cs, function(u) 
        sample(u, min(length(u), as.numeric(c)))))
}

# run method & write results to .rds
source(fun <- snakemake@input$fun)
fun <- gsub("(.R)", "", basename(fun))
res <- get(fun)(sim[gs, cs], meth_pars)

if (!inherits(res$tbl, "error")) {
    # add metadata
    gi <- metadata(sim)$gene_info %>% 
        dplyr::filter(gene %in% gs) %>% 
        dplyr::mutate_at("cluster_id", as.character) %>% 
        dplyr::select(-"logFC") %>% 
        dplyr::mutate(., sim_lfc = eval(parse(text = ifelse(
            .$sim_mean.B == 0, "0", "log2(sim_mean.B/sim_mean.A)"))))
    
    res$tbl <- left_join(gi, res$tbl, by = c("gene", "cluster_id")) %>%
        {if ("logFC" %in% names(.)) 
            dplyr::rename(., est_lfc = logFC) else .} %>%
        dplyr::mutate(did, sid, mid, i, j, g, c, k, s,
            is_de = as.integer(!category %in% c("ee", "ep")))
}

saveRDS(res, snakemake@output$res)