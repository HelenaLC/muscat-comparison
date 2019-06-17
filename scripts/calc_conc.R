# load packages
suppressMessages(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(purrr)
        library(reshape2)
    }))

# source utils
source(snakemake@input$utils)

#sim <- list.files("data/sim_data/kang", "de10;", full.names = TRUE)
#res <- list.files("results/kang/de10", full.names = TRUE)
gis <- .get_gis(snakemake@input$sim)
df <- .tidy_res(snakemake@input$res, gis)

# between-method concordance ---------------------------------------------------
ranks <- df %>% 
    mutate_at("p_val", function(u) { u[is.na(u)] <- 1; return(u) }) %>% 
    group_by(sim_rep, method, cluster_id) %>% 
    mutate(rank = order(p_val)) %>% 
    group_by(sim_rep, cluster_id) %>% 
    group_split %>% 
    map(acast, gene ~ method, value.var = "rank") 

mids <- names(method_colors)
mids <- expand.grid(mids, mids, stringsAsFactors = FALSE)
mids <- mids[apply(mids, 1, function(u) u[1] != u[2]), ]
ks <- vapply(gis, function(u) sum(u$is_de), numeric(1))
ks <- rep.int(ks, vapply(gis, function(u) 
    nlevels(factor(u$cluster_id)), numeric(1)))
cmt <- lapply(seq_along(ranks), function(i) {
    rs <- ranks[[i]]
    apply(mids, 1, function(ms) {
        if (!all(ms %in% colnames(rs)))
            return(NULL)
        ocs <- .calc_occs(rs[, ms], ks[i])
        auc <- .calc_aucs(ocs)
        auc$method1 <- ms[1]
        auc$method2 <- ms[2]
        return(auc)
    }) %>% bind_rows
})

mat <- lapply(seq_along(ranks), function(i) 
    filter(cmt[[i]], k == ks[i])) %>% 
    map(filter, nbr_occ == nbr_cols) %>% 
    map(dcast, method1 ~ method2, mean, value.var = "AUCs") %>%
    map(data.frame, row.names = 1, check.names = FALSE) %>% 
    map(function(u) { diag(u) <- 1; return(u) }) %>% 
    Reduce(f = "+") / length(ranks[[1]])

stopifnot(do.call("identical", dimnames(mat)))
stopifnot(all((mat == t(mat))[!is.na(mat == t(mat))]))

d <- 1 - mat
d[is.na(d)] <- 1
tree_avg <- hclust(as.dist(d))
subc_avg <- .get_subcl(tree_avg)
plot(tree_avg)

# calculate concordance ========================================================

# within-method concordance ----------------------------------------------------

# between-method concordance ---------------------------------------------------

# get subclusters for ea. instance ---------------------------------------------

subc_all <- lapply(seq_along(ranks), function(i) {
    mat <- acast(cmt[[i]], method1 ~ method2, mean, value.var = "AUCs") 
    diag(mat) <- 1 # add concordance for ea. method w/ itself
    stopifnot(identical(rownames(mat), colnames(mat)))
    stopifnot(all((mat == t(mat))[!is.na(mat == t(mat))]))
    d <- 1 - mat
    d[is.na(d)] <- 1
    h <- hclust(as.dist(d))
    .get_subcl(h)
})

# get stability values for ea. subcluster in subcl_avg
stab_sco <- rowMeans(sapply(subc_all, function(u) subc_avg %in% u))

res <- list(
    tree_avg = tree_avg, 
    subc_avg = subc_avg, 
    subc_all = subc_all, 
    stab_sco = stab_sco)

saveRDS(res, snakemake@output$res)