suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(purrr)
    library(SingleCellExperiment)
})

sim <- readRDS(args$sim)
meth_pars <- as.list(fromJSON(args$meth_pars))
run_pars <- as.list(fromJSON(args$run_pars))

set.seed(run_pars$seed + as.numeric(wcs$j))

# subset clusters & samples
kids <- levels(sim$cluster_id)
sids <- levels(sim$sample_id)
m <- match(sids, sim$sample_id)
gids <- sim$group_id[m]

if (wcs$k != "x") kids <- sample(kids, wcs$k)
if (wcs$s != "x") sids <- sapply(split(sids, gids), sample, wcs$s)

sim <- .filter_sce(sim, kids, sids)

# subset genes & cells
gs <- rownames(sim)
cs <- colnames(sim)

if (wcs$g != "x") 
    gs <- sample(gs, min(nrow(sim), as.numeric(wcs$g)))

if (wcs$c != "x") {
    cs <- split(cs, list(sim$cluster_id, sim$sample_id))
    cs <- unlist(lapply(cs, function(u) 
        sample(u, min(length(u), as.numeric(wcs$c)))))
}

# run method & write results to .rds
source(fun <- args$fun)
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
        dplyr::mutate(did = wcs$did, sid = wcs$sid, mid = wcs$mid, 
            i = wcs$i, j = wcs$j, g = wcs$g, c = wcs$c, k = wcs$k, s = wcs$s,
            is_de = as.integer(!category %in% c("ee", "ep")))
}

saveRDS(res, args$res)