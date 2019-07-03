config <- yaml::read_yaml("config.yaml")

ds10 <- c(0.9, 0, 0.1, 0, 0, 0)

sim_pars <- list(
    nill = list(nr = 3, p_dd = diag(6)[1, ], seed = 1),
    ds10 = list(nr = 5, p_dd = c(0.9, 0, 0.1, 0, 0, 0), seed = 10),
    dp10 = list(nr = 5, p_dd = c(0.9, 0, 0, 0.1, 0, 0), seed = 20),
    dm10 = list(nr = 5, p_dd = c(0.9, 0, 0, 0, 0.1, 0), seed = 30),
    db10 = list(nr = 5, p_dd = c(0.9, 0, 0, 0, 0, 0.1), seed = 40),
    
    ds10_nc = list(nk = 2, ns = 3, seed = 50, nc = 2*2*3*400),
    ds10_ns = list(nk = 2, ns = 5, seed = 60)
)

def_pars <- list(nr = 1, nk = 2, ns = 3, 
    ng = 2e3, nc = function(nk, ns) 2*nk*ns*100, 
    p_dd = ds10, seed = 1)

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
fns <- file.path(config$sim_pars, fns)
names(fns) <- sim_ids

for (id in sim_ids) {
    new <- sim_pars[[id]]
    if (file.exists(fns[id])) {
        old <- jsonlite::fromJSON(fns[id])
        if (!isTRUE(all.equal(old, new, tolerance = 1e-6)))
            write(jsonlite::toJSON(new, null = "null"), fns[id])
    } else {
        write(jsonlite::toJSON(new, null = "null"), fns[id])
    }
}
