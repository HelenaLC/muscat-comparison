args <- R.utils::commandArgs(
    trailingOnly = TRUE, 
    asValues = TRUE)

suppressMessages(library(muscat))
source(file.path("scripts", "utils.R"))

# load reference SCE
ref <- readRDS(args$input_sce)

# prep. SCE for simulation w/ 'muscat::simData'
sce <- prepSim(ref, verbose = FALSE,
    min_count = 1, min_cells = 10,   # keep genes w/ count > 1 in >= 10 cells
    min_genes = 100, min_size = 100) # keep cells w/ >= 100 detected genes

# write SCE to .rds
saveRDS(sce, args$output_sce)