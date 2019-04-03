# load packages
suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(jsonlite)
        library(magrittr)
    })
)

# specify simulation parameters
de10 <- c(.9,0,.1,0,0,0)
sim_ids <- names(
    sim_pars <- list(
        nill = list(n_reps = 3, n_genes = 500, n_cells = 12*200, p_dd = diag(6)[1, ], seed = 1, nk = 2, ns = 3),
        
        de10 = list(n_reps = 3, n_genes = 500, n_cells = 12*200, p_dd = de10,             seed = 10, nk = 2, ns = 3),
        dp10 = list(n_reps = 3, n_genes = 500, n_cells = 12*200, p_dd = c(.9,0,0,.1,0,0), seed = 20, nk = 2, ns = 3),
        dm10 = list(n_reps = 3, n_genes = 500, n_cells = 12*200, p_dd = c(.9,0,0,0,.1,0), seed = 30, nk = 2, ns = 3),
        db10 = list(n_reps = 3, n_genes = 500, n_cells = 12*200, p_dd = c(.9,0,0,0,0,.1), seed = 40, nk = 2, ns = 3),
        
        #nill_ng = list(n_reps = 1, n_genes = 3e3, n_cells = 12*200, p_dd = diag(6)[1, ], seed = 5, nk = 2, ns = 3),
        de10_nc = list(n_reps = 1, n_genes = 500, n_cells = 12*250, p_dd = de10, seed = 50, nk = 2, ns = 3),
        de10_ns = list(n_reps = 1, n_genes = 500, n_cells = 16*200, p_dd = de10, seed = 60, nk = 2, ns = 4),
        
        mixed_ss = list(n_reps = 3, n_genes = 500, n_cells = 12 * 100, p_dd = de10, probs = list(NULL, c(1/9,3/9,5/9), NULL), seed = 70, nk = 2, ns = 3)
        #unbal_gs = list(n_genes = 500, n_cells = 6 * 250, p_dd = de10, fc = 1.5, probs = list(NULL, NULL, c(.3,.7)),    seed = 7, nk = 2, ns = 3)
    )
)

# create result directories
dir_res <- file.path("results", with(
    expand.grid(snakemake@config$data_ids, sim_ids),
    file.path(Var1, Var2)))
for (dir in dir_res)
    if (!dir.exists(dir)) dir.create(dir)

# write simulation IDs to .txt
write(sim_ids, snakemake@config$sim_ids)

# write to .json (only if something changed!)
fns <- sprintf("%s.json", sim_ids)
fns <- file.path(snakemake@config$sim_pars, fns) %>% set_names(sim_ids)
for (id in sim_ids) {
    if (file.exists(fns[id])) {
        pars <- fromJSON(fns[id])
        if (!isTRUE(all.equal(pars, sim_pars[[id]])))
            write(toJSON(sim_pars[[id]], null = "null"), fns[id])
    } else {
        write(toJSON(sim_pars[[id]], null = "null"), fns[id])
    }
}

# run-mode configurations ------------------------------------------------------
L <- list(
    nill_ng = list(seed = 31, n_reps = 5, n_genes = c(100,250,500,1e3,2e3)),
    de10_nc = list(seed = 67, n_reps = 3, n_cells = c(10,25,50,100,200)),
    de10_ns = list(seed = 92, n_reps = 3, n_samples = c(2, 3, 4))
)

# default config.
default_run_pars <- list(seed = 1, n_reps = 1, n_genes = "all", n_cells = "all", n_samples = "all")
L <- lapply(L, function(l) {
    missing <- !names(default_run_pars) %in% names(l)
    l[names(default_run_pars)[missing]] <- default_run_pars[missing]
    return(l)
})

missing <- sim_ids[!sim_ids %in% names(L)]
L[missing] <- replicate(length(missing), default_run_pars, simplify = FALSE)

# write to .json (only if something changed!)
fns <- paste0(sim_ids, ".json")
fns <- file.path("config", fns) %>% set_names(sim_ids)
for (id in sim_ids) {
    if (file.exists(fns[id])) {
        l <- fromJSON(fns[id])
        if (!isTRUE(all.equal(l, L[[id]])))
            write(toJSON(L[[id]]), fns[id])
    } else {
        write(toJSON(L[[id]]), fns[id])
    }
}
