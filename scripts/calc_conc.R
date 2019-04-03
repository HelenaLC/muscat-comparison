# load packages
suppressMessages(
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(purrr)
    }))

sim_dir <- file.path("data", "sim_data", "kang")
res_dir <- file.path("results", "kang")
sim_ids <- readLines("metadata/sim_ids.txt")[-1] %>% purrr::set_names()
sim_ids <- grep("de10", sim_ids, value = TRUE)
method_pars <- read.csv("metadata/method_ids.csv", row.names = 1)
method_ids <- as.character(method_pars$id)

# get rankings for each sim_id, sim_rep, run_rep, method
p_val <- lapply(sim_ids, function(sim_id) {
    sim <- list.files(sim_dir, paste0(sim_id, ";"), full.names = TRUE)
    res <- list.files(file.path(res_dir, sim_id), full.names = TRUE)
    
    # load truth & tidy results
    gis <- .get_gis(sim)
    .tidy_res(res, gis)
})
df <- p_val %>% bind_rows(.id = "sim_id") %>% 
    mutate_at("p_val", function(u) { u[is.na(u)] <- 1; return(u) }) %>% 
    group_by(sim_id, sim_rep, n_genes, n_cells, n_samples, run_rep, method) %>% 
    mutate(rank = order(p_val)) %>% ungroup %>% 
    data.frame(stringsAsFactors = FALSE)

# calculate concordance ========================================================
res <- list()

# within-method concordance ----------------------------------------------------
tmp <- df %>% setDT %>% split(
    flatten = FALSE, sorted = TRUE, by = c("sim_id", "sim_rep", 
        "n_genes", "n_cells", "n_samples", "method", "run_rep")) %>% 
    map_depth(-2, pull, "rank") %>% 
    map_depth(-2, function(u) {
        len <- vapply(u, length, numeric(1))
        data.frame(u[len != 0], 
            check.names = FALSE,
            stringsAsFactors = FALSE)
    }) %>% map_depth(-2, as.matrix)
while (!any(sapply(tmp, is, "matrix"))) 
    tmp <- flatten(tmp)
n_reps <- vapply(tmp, ncol, numeric(1))
tmp <- tmp[n_reps > 1]

auc_in <- lapply(tmp, function(u) {
    ocs <- calculate_nbr_occurrences(u, 100)
    calc_auc(ocs)
}) %>% bind_rows(.id = "id")

stab <- auc_in %>% 
    filter(k == 100) %>% 
    group_by(id) %>% 
    summarize(meanAUC = mean(AUCs))

# between-method concordance ---------------------------------------------------
tmp <- df %>% setDT %>% split(
    flatten = FALSE, sorted = TRUE, by = c("sim_id", "sim_rep", 
        "n_genes", "n_cells", "n_samples", "method")) %>% 
    map_depth(-2, pull, "rank") %>% 
    map_depth(-2, function(u) {
        len <- vapply(u, length, numeric(1))
        data.frame(u[len != 0], 
            check.names = FALSE,
            stringsAsFactors = FALSE)
    }) %>% map_depth(-2, as.matrix)
while (!any(sapply(tmp, is, "matrix"))) tmp <- flatten(tmp)
   
ms <- expand.grid(method_ids, method_ids, stringsAsFactors = FALSE) 
ms <- ms[apply(ms, 1, function(u) u[1] != u[2]), ]
auc_bw <- apply(ms, 1, function(ms) {
    lapply(tmp, function(u) {
        if (!all(ms %in% colnames(u))) return(NULL)
        ocs <- calculate_nbr_occurrences(u[, ms], 100)
        auc <- calc_auc(ocs)
        auc$method1 <- ms[1]
        auc$method2 <- ms[2]
        return(auc)
    })
}) %>% flatten %>% bind_rows

mat <- auc_bw %>% 
    filter(k == 100) %>% 
    group_by(method1, method2) %>% 
    summarize(meanAUC = mean(AUCs)) %>% 
    dcast(method1 ~ method2, value.var = "meanAUC") %>% 
    data.frame(row.names = 1, check.names = FALSE)
#mat[is.na(mat)] <- 1

# hierarchical clustering ------------------------------------------------------
d <- 1 - mat
d[is.na(d)] <- 1
d <- as.dist(d)
h <- hclust(d)
plot(h)

ggt <- ggtree(tidytree::as.phylo(h))
subcl <- get_subclusters(h)

# for (m in seq_len(length(subcl))) {
#     tryCatch({
#         i <- MRCA(ggt, subcl[[m]])
#         if (stab[m] >= 0.1)
#             ggt$data[i, "label"] <- round(stab[m], 2)
#     }, error = function(e) NULL)
# }

ggt <- ggt + geom_label2(aes(subset = !isTip, label = label), size = 2) +
    geom_tiplab(aes(angle = 90), hjust = 1) +
    ggplot2::scale_x_reverse() + ggplot2::coord_flip() +
    theme(plot.margin = unit(c(0, 0, 10, 0), "mm")) +
    xlim_tree(1)

ss <- gsub("\\(.*", "", method_ids)
ss[method_pars$type != "pb"] <- NA

inp <- gsub("\\..*", "", method_pars$id) %>% set_names(method_pars$id)
inp[grep("mm|scDD|AD", inp)] <- "logcounts"
inp[grep("MAST", inp)] <- "logcpm"
typ <- method_pars$type
levels(typ) <- paste0(c("group", "group", "cell", "sample", "group"), "-level")
anno <- method_pars %>% 
    mutate(
        type = typ,
        method = id,
        log = as.logical(c(0,0,0,0,0,0,0,1,0,0,1,0,1,1,1,1,1,1,1,1)),
        input = factor(inp),
        aggregation = ss)

tiporder <- ggt$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
anno <- anno[match(tiporder$label, anno$method), ]

cols <- list(
    type = setNames(brewer.pal(3, "Reds"),
        paste0(c("cell", "sample", "group"), "-level")),
    input = setNames(colorRampPalette(brewer.pal(9, "Blues"))(nlevels(anno$input)), levels(anno$input)),
    log = setNames(c("grey80", "grey60"), c("TRUE", "FALSE")),
    aggregation = setNames(
        brewer.pal(3, "Oranges"),
        c("sum", "mean", "median")))

hms <- lapply(seq_along(cols), function(i) {
    ggplot(anno) + geom_tile(aes_string(x = seq_len(nrow(anno)), y = 1, fill = names(cols)[i])) + 
        geom_vline(xintercept = (seq_len(nrow(anno)-1)) + 0.5, 
            linetype = "solid", color = "white", size = 0.5) + 
        scale_fill_manual(values = cols[[i]]) + 
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + 
        theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(0, 5, 0, 5), "mm"))
})

ggt <- ggt + theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
lgd <- lapply(hms, get_legend)
lgd <- plot_grid(plotlist = lgd, ncol = 1)
hms <- lapply(hms, function(p) p + theme(legend.position = "none"))
plot_grid(
    plot_grid(
        plotlist = c(list(ggt), hms),
        ncol = 1, rel_heights = c(5, rep(1, length(cols)))),
    lgd, ncol = 2, 
    rel_widths = c(1, 0.5))

## For each k in 1:maxrank, count the number of genes occurring each number of 
## times. The output is a data frame with three columns: nbr_occ, nbr_genes, k. 
## For a given row, interpret as follows: among the top-k genes from each 
## method, nbr_genes occur exactly nbr_occ times. 
calculate_nbr_occurrences <- function(mtx, maxrank) {
    maxrank <- min(maxrank, nrow(mtx))
    if (ncol(mtx) > 1) {
        M <- matrix(0, max(mtx[1:maxrank, ]), maxrank)
        for (i in 1:ncol(mtx)) {
            M[cbind(mtx[1:maxrank, i], 1:maxrank)] <- M[cbind(mtx[1:maxrank, i], 1:maxrank)] + 1
        }
        M <- M[rowSums(M) != 0, ]
        M <- t(apply(M, 1, cumsum))
        M2 <- matrix(0, ncol(mtx), maxrank)
        for (i in 1:nrow(M)) {
            M2[cbind(M[i, ], 1:ncol(M))] <- M2[cbind(M[i, ], 1:ncol(M))] + 1
        }
        M2 <- M2 %>% reshape2::melt() %>%
            dplyr::rename(k = Var2, nbr_genes = value, nbr_occ = Var1)
        M2 %>% dplyr::arrange(k, nbr_occ) %>%
            dplyr::mutate(nbr_cols = ncol(mtx))
    } else {
        NULL
    }
}

## Calculate partial (cumulative) AUCs. 
## Assumes that x variable = k, y variable = nbr_genes
calc_auc <- function(x) {
    x %>% dplyr::mutate(dx = c(k[1], diff(k)),
        dy = c(nbr_genes[1], diff(nbr_genes)),
        ys = c(0, nbr_genes[-length(nbr_genes)])) %>%
        dplyr::mutate(AUC = cumsum(dx * dy/2 + dx * ys)) %>%
        dplyr::mutate(AUCs = AUC/(k^2/2))
}

## Get all subclusters from an hclust object
get_subclusters <- function(hcl) {
    m <- hcl$merge
    labs <- hcl$labels
    L <- list()
    for (i in seq_len(nrow(m))) {
        tmp <- c()
        if (m[i, 1] < 0) tmp <- c(tmp, labs[-m[i, 1]])
        else tmp <- c(tmp, L[[m[i, 1]]])
        if (m[i, 2] < 0) tmp <- c(tmp, labs[-m[i, 2]])
        else tmp <- c(tmp, L[[m[i, 2]]])
        L[[i]] <- sort(tmp)
    }
    L
}
