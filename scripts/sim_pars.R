config <- yaml::read_yaml("config.yaml")

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
        nill = list(n_reps = 3, n_genes = 5e3, n_cells = 12*200, p_dd = diag(6)[1, ], seed = 1),
        
        de10 = list(n_reps = 10, n_genes = 4e3, n_cells = 12*200, p_dd = de10,             seed = 10),
        dp10 = list(n_reps = 10, n_genes = 4e3, n_cells = 12*200, p_dd = c(.9,0,0,.1,0,0), seed = 30),
        dm10 = list(n_reps = 10, n_genes = 4e3, n_cells = 12*200, p_dd = c(.9,0,0,0,.1,0), seed = 40),
        db10 = list(n_reps = 10, n_genes = 4e3, n_cells = 12*200, p_dd = c(.9,0,0,0,0,.1), seed = 60),
        
        #nill_ng = list(n_reps = 1, n_genes = 3e3, n_cells = 12*200, p_dd = diag(6)[1, ], seed = 5),
        #de10_nc = list(n_reps = 10,  n_genes = 5e3, n_cells =  8*500, p_dd = de10, seed = 50, nk = 2, ns = 2),
        de10_ns = list(n_reps = 10, n_genes = 4e3, n_cells = 20*200, p_dd = de10, ns = 5, seed = 80),
        
        mixed_ss = list(n_reps = 10, n_genes = 4e3, n_cells = 12*200, p_dd = de10, probs = list(NULL, c(1/9,3/9,5/9), NULL), seed = 100)
        #unbal_gs = list(n_genes = 500, n_cells = 6 * 250, p_dd = de10, fc = 1.5, probs = list(NULL, NULL, c(.3,.7)),    seed = 7)
    )
)

# create result directories
dir_res <- file.path("results", with(
    expand.grid(config$dat_ids, sim_ids),
    file.path(Var1, Var2)))
for (dir in dir_res)
    if (!dir.exists(dir)) dir.create(dir)

# write simulation IDs to .txt
write(sim_ids, config$sim_ids)

# write to .json (only if something changed!)
fns <- sprintf("%s.json", sim_ids)
fns <- file.path(config$sim_pars, fns) %>% set_names(sim_ids)
for (id in sim_ids) {
    if (file.exists(fns[id])) {
        old_pars <- fromJSON(fns[id])
        if (!isTRUE(all.equal(old_pars, sim_pars[[id]], tolerance = 1e-4)))
            write(toJSON(sim_pars[[id]], null = "null"), fns[id])
    } else {
        write(toJSON(sim_pars[[id]], null = "null"), fns[id])
    }
}

# run-mode configurations ------------------------------------------------------)
L <- list(
    kang = list(de10_ns = list(seed = 92, n_reps = 1, n_samples = c(2, 3, 4, 5))), 
    magl = list(de10_ns = "norun")
)

# default config.
default_run_pars <- list(seed = 1, n_reps = 1, n_genes = "all", n_cells = "all", n_samples = "all")
L <- lapply(L, lapply, function(l) {
    if (l[1] == "norun") return(l)
    missing <- !names(default_run_pars) %in% names(l)
    l[names(default_run_pars)[missing]] <- default_run_pars[missing]
    return(l)
})
L <- lapply(L, function(l) {
    missing <- sim_ids[!sim_ids %in% names(l)]
    l[missing] <- replicate(length(missing), default_run_pars, simplify = FALSE)
    return(l)
})

# write to .json (only if something changed!)
fns0 <- paste0(sim_ids, ".json")
for (dat_id in names(L)) {
    fns <- paste(dat_id, fns0, sep = ";")
    dirs <- file.path("config", fns) %>% set_names(sim_ids)
    for (sim_id in sim_ids) {
        if (file.exists(dirs[sim_id])) {
            l <- fromJSON(dirs[sim_id])
            if (!isTRUE(all.equal(l, L[[dat_id]][[sim_id]])))
                write(toJSON(L[[dat_id]][[sim_id]]), dirs[sim_id])
        } else {
            write(toJSON(L[[dat_id]][[sim_id]]), dirs[sim_id])
        }
    }
}
