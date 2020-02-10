config <- yaml::read_yaml("config.yaml")

de10 <- c(0.9, 0, 0.1, 0, 0, 0)

sim_pars <- list(
    # null simulations
    nill = list(nr = 2, p_dd = diag(6)[1, ], seed = 1), 
    # pure simulations
    de10 = list(nr = 2, p_dd = c(0.9, 0, 0.1, 0, 0, 0), seed = 10),
    dp10 = list(nr = 2, p_dd = c(0.9, 0, 0, 0.1, 0, 0), seed = 30),
    dm10 = list(nr = 2, p_dd = c(0.9, 0, 0, 0, 0.1, 0), seed = 50),
    db10 = list(nr = 2, p_dd = c(0.9, 0, 0, 0, 0, 0.1), seed = 70),
    # increasing number of genes, cells, replicates
    de10_ng = list(nr = 2, nk = 2, ns = 3, seed = 80, nc = 2*2*3*100),
    de10_nc = list(nr = 2, nk = 2, ns = 3, seed = 90, nc = 2*2*3*500),
    de10_ns = list(nr = 2, nk = 2, ns = 5, seed = 110),
    # mixed
    dx = list(nr = 2, p_dd = c(0.8, 0, rep(0.05, 4)), seed = 200)
)

# unbalanced sample sizes ------------------------------------------------------
ss_ns <- 3
ss <- lapply(seq_len(4), function(i) {
    ss <- rep(1, ss_ns) / seq(1, i, length = ss_ns)
    ss / sum(ss)
})
for (i in seq_along(ss)) {
    id <- paste0("de10_ss", i)
    sim_pars[[id]] <- list(
        nr = 2, nk = 2, ns = ss_ns, seed = 130+20*(i-1),
        p_dd = de10, probs = list(NULL, ss[[i]], NULL))
}
# unbalanced group sizes -------------------------------------------------------
gs <- list(c(0.5, 0.5), c(0.45, 0.55), c(0.4, 0.6), c(0.3, 0.7))
for (i in seq_along(gs)) {
    id <- paste0("de10_gs", i)
    sim_pars[[id]] <- list(
        nr = 2, nk = 2, seed = 150+20*(i-1),
        p_dd = de10, probs = list(NULL, NULL, gs[[i]]))
}

def_pars <- list(nr = 1, nk = 3, ns = 3, 
    ng = 2e3, nc = function(nk, ns) 2*nk*ns*200, 
    p_dd = de10, probs = NULL, seed = 1)

sim_pars <- lapply(sim_pars, function(u) {
    v <- !names(def_pars) %in% names(u)
    u[names(def_pars)[v]] <- def_pars[v]
    if (is.function(u$nc))
        u$nc <- u$nc(u$nk, u$ns)
    return(u)
})

sim_ids <- names(sim_pars)
write(jsonlite::toJSON(sim_ids), config$sids)

# write parameters to .json (only if something changed!)

fns <- paste0(sim_ids, ".json")
fns <- paste0(config$sim_pars, fns)
names(fns) <- sim_ids

for (id in sim_ids) {
    new <- sim_pars[[id]]
    if (file.exists(fns[id])) {
        old <- jsonlite::fromJSON(fns[id])
        if (!isTRUE(all.equal(old, new, tolerance = 1e-3)))
            write(jsonlite::toJSON(new, null = "null"), fns[id])
    } else {
        write(jsonlite::toJSON(new, null = "null"), fns[id])
    }
}
