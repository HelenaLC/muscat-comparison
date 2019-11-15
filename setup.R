config <- yaml::read_yaml("config.yaml")
x <- c("sim_pars", "run_pars", "meth_pars", "raw_data", 
    "sim_data", "results", "plots", "figures", "logs")
for (dir in unlist(config[x]))
    if (!dir.exists(dir) & !isTRUE(grep("\\.", dir) == 1))
        dir.create(dir, recursive = TRUE)

for (script in file.path(config$scripts, 
    paste0(c("sim_pars", "run_pars", "meth_pars"), ".R")))
    source(script)