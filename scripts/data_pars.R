config <- yaml::read_yaml("config.yaml")

# for each dataset, create a list of
#   dir_raw:        directory to raw counts
#   dir_gene_md:    directory to gene metadata (SCE rowData)
#   dir_cell_md:    directory to cell metadata (SCE colData)
#   filter_cells:   NULL or function to filter cells from colData 
#   sample_id, group_id, cluster_id:
#       the colData column containing xx IDs or
#       a function to obtain xx IDs from the colData
#   group_keep:     group to keep for simulation

data_ids <- names(
    data_pars <- list(
        kang = list(
            dir_raw = c("GSM2560248_2.1.mtx", "GSM2560249_2.2.mtx"),
            dir_gene_md = "GSE96583_batch2_genes.tsv",
            dir_cell_md = "GSE96583_batch2_df.tsv",
            filter_cells = function(x) x$multiplets == "singlet" & !is.na(x$cell),
            sample_id = function(x) factor(with(x, paste0(stim, ind))),
            group_id = "stim",
            cluster_id = "cell",
            group_keep = "ctrl"
        ),
        magl = list(
            dir_raw = "MAGL-SCE.rds",
            dir_gene_md = NULL,
            dir_cell_md = NULL,
            sample_id = "sample_id",
            group_id = "group_id",
            cluster_id = "cluster_id",
            group_keep = "WT"
        )
    )
)

# prefix directories
data_pars <- setNames(lapply(data_ids, function(id) {
    dirs <- grep("dir", names(data_pars[[id]]))
    data_pars[[id]][dirs] <- lapply(data_pars[[id]][dirs], 
        function(u) file.path(config$raw_data, id, u))
    return(data_pars[[id]])
}), data_ids)

# create directories
for (u in c("sim_data", "results", "figures")) {
    base <- config[[u]]
    dirs <- file.path(base, data_ids)
    for (dir in dirs) 
        if (!dir.exists(dir)) 
            dir.create(dir)
}

# write parameters to .rds
# (only if something changed!)
dir <- config$data_pars
if (!dir.exists(dir))
    dir.create(dir)

fns <- sprintf("%s.rds", file.path(dir, data_ids))
names(fns) <- data_ids
for (id in data_ids) {
    fn <- fns[[id]]
    if (file.exists(fn)) {
        pars <- readRDS(fn)
        if (!identical(pars, data_pars[[id]]))
            saveRDS(data_pars[[id]], fn)
    } else {
        saveRDS(data_pars[[id]], fn)
    }
}
