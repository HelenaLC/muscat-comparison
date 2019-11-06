config <- yaml::read_yaml("config.yaml")
x <- c("sim_pars", "run_pars", "meth_pars", 
    "raw_data" "sim_data", "results", "figures")
for (dir in unlist(config[x]))
    if (!dir.exists(dir) & !isTRUE(grep("\\.", dir) == 1))
        dir.create(dir, recursive = TRUE)

scripts <- file.path(config$scripts,
    paste0(c("sim_pars", "run_pars", "meth_pars"), ".R"))
sapply(scripts, source)