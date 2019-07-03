config <- yaml::read_yaml("config.yaml")
sim_ids <- jsonlite::fromJSON(config$sids)

def_pars <- list(nr = 1, nk = "x", ns = "x", ng = "x", nc = "x", seed = 1)

run_pars <- list(
    kang = list(
        ds10_nc = list(nr = 3, nc = c(50, 100, 200, 300)),
        ds10_ns = list(nr = 3, ns = c(2, 3, 4, 5))
    ),
    magl = list(
        ds10_nc = list(nr = 3, nc = c(50, 100, 200, 300)),
        ds10_ns = NULL
    )
)

run_pars <- lapply(run_pars, lapply, function(l) {
    if (is.null(l)) return(l)
    i <- !names(def_pars) %in% names(l)
    l[names(def_pars)[i]] <- def_pars[i]
    return(l)
})

run_pars <- lapply(run_pars, function(l) {
    i <- sim_ids[!sim_ids %in% names(l)]
    l[i] <- replicate(length(i), def_pars, simplify = FALSE)
    return(l)
})

# write parameters to .json (only if something changed!)

fns0 <- paste0(sim_ids, ".json")
for (data_id in config$dids) {
    fns <- paste(data_id, fns0, sep = ";")
    fns <- file.path(config$run_pars, fns)
    names(fns) <- sim_ids
    for (sim_id in sim_ids) {
        new <- run_pars[[data_id]][[sim_id]]
        if (file.exists(fns[sim_id])) {
            old <- jsonlite::fromJSON(fns[sim_id])
            if (!isTRUE(all.equal(old, new)))
                write(jsonlite::toJSON(new, null = "null"), fns[sim_id])
        } else {
            write(jsonlite::toJSON(new, null = "null"), fns[sim_id])
        }
    }
}
