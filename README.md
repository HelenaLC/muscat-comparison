# On the discovery of population-specific state transitions <br> from multi-sample multi-condition scRNA-seq data

This repository contains all the necessary code to perform the evaluations and analyses from our preprint available on [bioRxiv](https://www.biorxiv.org/content/10.1101/713412v1).

Analyses discussed in the **Differential state analysis of mouse cortex exposed to LPS treatment** results section are provided as a browsable `workflowr`<sup>[1](#f1)</sup> website [HERE](http://htmlpreview.github.io/?https://github.com/HelenaLC/muscat-comparison/blob/master/LPS/docs/index.html).

## Prerequisites

For installation of the required libraries, we'll fist install the `r BiocStyle::Biocpkg("BiocManager")` package:

```r
install.packages("BiocManager")
```

The code in this repository was developed using **R v3.6.2** and **Bioconductor v3.10**. Versions of R and Bioconductor that are currently being run should be checked via:

```r
version
BiocManager::version()
```

Finally, the code chunk below will install all package dependencies:

```r
# install 'ggrastr' from GitHub
BiocManager::install("VPetukhov/ggrastr")

# install packages from CRAN & Bioconductor
pkgs <- c("AnnotationDbi","circlize","countsimQC","cowplot","data.table",
    "DESeq2","DropletUtils","dplyr","edgeR","ggplot2","iCOBRA","kSamples",
    "jsonlite","limma","M3C","magrittr","MAST","Matrix","muscat","msigdbr",
    "org.Mm.eg.db","pheatmap","purrr","RColorBrewer","readxl","reshape2",
    "S4Vectors","scater","scDD","scds","scran","sctransform","Seurat",
    "SingleCellExperiment","topGO","UpSetR","viridis","workflowr","yaml")
BiocManager::install(pkgs, ask = FALSE)
```

## Setup

R version and library have to be specified under `R` in the `config.yaml` file (e.g., `R: "R_LIBS_USER=/path/to/library /path/to/R/executable"`). If you run into any issues, I recommend running the specified character string from the command line and assuring that the outputs of `version` and `.libPaths()` are what you expect them to be.

Without modifications, the `Snakemake` comparison relies on 2 reference datasets for data simulation, method execution and comparison. These (or any other references) have to be downloaded, saved as `.rds` objects with appropriate names, and placed inside the `data/raw_data` directory. This is exemplified here for the Kang et al. reference used in the preprint:

```r
# install & load 'ExperimentHub'
BiocManager::install("ExperimentHub")
library(ExperimentHub)

# initialize hub instance
eh <- ExperimentHub()

# list data available in 'muscData'
(q <- query(eh, c("Kang", "muscData")))

# load 'SingleCellExperiment's using IDs from above
sce <- eh[[q$ah_id]]

# save as .rds; name should be '<id>_sce0.rds'
fn <- "kang_sce0.rds"
dir <- file.path("...", "data", "raw_data")
saveRDS(sce, file.path(dir, fn))
```

Finally, execution of the `Snakemake` file requires running the `setup.R` script **once** to create all required directories as well as simulation, method and run parameters:

```r
# from within R
source("setup.R")

# from the terminal
Rscript setup.R
```

The `Snakemake` should run now. A couple more points to note:

1. `sim/run/meth_pars.R` in the `scripts` are re-exected with every `Snakemake` run, and any changes made to them will automatically be recognized (e.g., when a new simulation scenario or method is added).
1. Running the whole workflow is computationally expensive (~3 days using 40 cores). For development purposes, with recommend limiting to 1 reference, fewer simulation replicates and/or fewer genes per simulation. Most importantly, at least initially, including one or no mixed model based methods will greatly speed things up!

## How to...

1. **add a new reference**
    * `<id>_sce0.rds` has to be in place as described above
    * `"<id>"` has to be added under `dids` in the `config.yaml` file
    * a `scripts/prep_<id>.R` has to be added to, for example, assure unique sample identifiers exist, remove un-assigned cells or cell multiplets etc.
1. **skip an existing reference**
    * simply remove the corresponding ID under `dids` in `config.yaml`
1. **add a new method**
    * for a single method, add a new line (with unique identifier) for that method in the corresponding `data.frame` constructed in `scripts/meth_pars.R`
    * for a new group of methods, add `id` under `ids` in the first line of `scripts/meth_pars.R` and code to construct a `data.frame` of appropriate format (must include a `id` column; see current methods for examples). Secondly, add a `apply_<id>.R` script under `scripts` that takes as input a SCE with `colData` columns `cluster/sample/group_id` and returns a `data.frame` with `p_adj.loc` and `p_adj.glb` values for each cluster-gene (see current `apply_x.R` scripts for exmples)
1. **skip an existing method**
    * to exclude a group of methods, comment out the corresponding `ids` in the first line of `scripts/meth_pars.R` (e.g., to skip all mixed model based methods, one would comment out `"mm"`)
    * to exclude a single method, remove that method from the corresponding `data.frame` in `scripts/meth_pars.R`

***

### Workflow structure in detail

In brief, our `Snakemake` workflow for method comparison is organized into

- a `config.yaml` file specify key parameters and directories
- a `scripts` folder housing all utilized scripts (see below)
- a `data` folder containing raw (reference) and simulated data
- a `meta` folder for simulation, runmode, and method parameters
- a `results` folder where all results are generated (as `.rds` files)
- a `plots` folder where all output plots are generated  
(as `.pdf` or `.png` and `.rds` files for `ggplot` objects)

The table below summarizes the different R scripts in `scripts`:

script      | description 
:-----------|:-----------------------------------------------
`prep_X`    | generates a references SCE for simulation by<br>i) keeping samples from one condition only; and,<br>ii) unifying relevant cell metadata names to `"cluster/sample/group_id"`
`prep_sim` | prepares a reference SCE for simulation by<br>i) retaining subpopulation-sample combinations with at least 100 cells; and,<br>ii) estimating cell / gene parameters (offsets / coefficients and dispersions)
`sim_pars`  | for ea. simulation ID, generates a `.json` file in `meta/sim_pars`<br>that specifies simulation parameters (e.g., prob. of DS, nb. of simulation replicates)
`run_pars`  | for ea. reference and simulation ID, generates a `.json` file in `meta/run_pars`<br>that specifies runmode parameters (e.g., nb. of cells/genes to sample, nb. of run replicates) 
`meth_pars` | for ea. method ID, generates a `.json` file in `meta/meth_pars`<br>that specifies method parameters
`sim_data`  | provided with a reference dataset and simulation parameters,<br>simulates data and writes a SCE to `data/sim_data`
`apply_X`   | wrapper to run DS method of type X (`pb`, `mm`, `ad`, `mast`, `scdd`)
`run_meth`  | reads in simulated data, method parameters, and performs DS analysis<br>by running the corresponding `apply_X` script
`run_meth_lps` | wrapper to apply method to the LPS dataset
`plot_null` | for ea. reference ID, plots nominal p-value distributions for all null simulations
`plot_perf_cat`     | plots TPR-FDR-points across DD categories for ea. p-value adjustment type (`p_adj.loc/glb`)
`plot_perf_by_nx`   | plots TPR-FDR-points across the nb. of `x` (cells = `c`, samples = `s`)
`plot_perf_by_xs`   | plots TPR-FDR-points across increasingly unbalanced sample/group-sizes
`plot_perf_by_expr` | plots TPR-FDR-points across expression-level groups
`plot_upset`        | plots an upset plot for the top gene-subpopulation combinations across methods and simulation replications
`plot_lfc`          | scatter plots of simulated vs. estimated logFC stratified by method and DD category
`plot_pb_mean_disp` | provided with a reference dataset, simulates a null dataset (no DS, no type-genes)<br>and plots pseudobulk-level mean-dispersion estimates for simulated vs. reference data
`plot_runtimes`     | barplots of runtimes vs. nb. of genes/cells
`utils`        | various helpers for data handling, formatting, and plotting
`session_info` | generates a `.txt` file capturing the output of `session_info()`

### References

<a name="f1">[1]</a>:
John Blischak, Peter Carbonetto and Matthew Stephens (2019).  
workflowr: A Framework for Reproducible and Collaborative Data Science.  
R package version 1.4.0. https://CRAN.R-project.org/package=workflowr