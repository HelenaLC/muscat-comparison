### On the discovery of population-specific state transitions <br> from multi-sample multi-condition scRNA-seq data

This repository contains all the necessary code to perform the evaluations and analyses from ...

### LPS dataset analysis

Analyses discussed in the **Differential state analysis of mouse cortex exposed to LPS treatment**  
results section are provided as a browsable `workflowr`<sup>[1](#f1)</sup> website [HERE](http://htmlpreview.github.io/?https://github.com/HelenaLC/muscat-comparison/blob/master/MAGL/docs/index.html).

### Comparison of DS analysis methods

In brief, our `snakemake` workflow for method comparison is organized into

- a `config.yaml` file specify key parameters and directories
- a `scripts` folder housing all utilized scripts (see below)
- a `data` folder containing raw (reference) and simulated data
- a `meta` folder for simulation, runmode, and method parameters
- a `results` folder where all results are generated (as `.rds` files)
- a `figures` folder where all output plots are generated  
(as `.pdf` or `.png` files, or `.rds` files for `ggplot` objects)

The table below summarizes the different R scripts in `scripts`:

script      | description 
------------|------------------------------------------------
`prep_X`    | generates a references SCE for simulation by<br>i)by keeping samples from one condition only; and,<br>ii) unifying relevant cell metadata names to `"cluster/sample/group_id"`
`prep_sim` | prepares a reference SCE for simulation by<br>i) retaining subpopulation-sample combinations with at least 100 cells; and,<br>ii) estimating cell / gene parameters (offsets / coefficients and dispersions)
`sim_pars`  | for ea. simulation ID, generates a `.json` file in `meta/sim_pars`<br>that specifies simulation parameters (e.g., prob. of DS, nb. of simulation replicates)
`run_pars`  | for ea. reference and simulation ID, generates a `.json` file in `meta/run_pars`<br>that specifies runmode parameters (e.g., nb. of cells/genes to sample, nb. of run replicates) 
`meth_pars` | for ea. method ID, generates a `.json` file in `meta/meth_pars`<br>that specifies method parameters
`sim_data`  | provided with a reference dataset and simulation parameters,<br>simulates data and writes a SCE to `data/sim_data`
`apply_X`   | wrapper to run DS method of type X (`pb`, `mm`, `ad`, `mast`, `scdd`)
`run_meth`  | reads in simulated data, method parameters, and performs DS analysis<br>by running the corresponding `apply_X` script
`plot_null` | for ea. reference ID, plots nominal p-value distributions for all null simulations
`plot_tprfdr`       | plots TPR-FDR-curves for a single result
`plot_perf_cat`     | plots TPR-FDR-points across DD categories for ea. p-value adjustment type (`p_adj.loc/glb`)
`plot_perf_by_nx`   | plots TPR-FDR-points across the nb. of `x` (cells = `c`, samples = `s`)
`plot_perf_by_ss`   | plots TPR-FDR-points across increasingly unbalanced sample-sizes
`plot_perf_by_expr` | plots TPR-FDR-points across expression-level groups
`plot_upset`        | plots an upset plot for the top gene-subpopulation combinations across methods and simulation replications
`plot_lfc`          | scatter plots of simulated vs. estimated logFC stratified by method and DD category
`plot_pb_mean_disp` | provided with a reference dataset, simulates a null dataset (no DS, no type-genes)<br>and plots pseudobulk-level mean-dispersion estimates for simulated vs. reference data
`plot_runtimes`     | barplots of runtimes vs. nb. of genes/cells
`utils`     | various helpers for data handling, formatting, and plotting

### References

<a name="f1">1</a>: John Blischak, Peter Carbonetto and Matthew Stephens (2019). workflowr: A Framework for Reproducible and Collaborative Data Science. R package version 1.4.0. https://CRAN.R-project.org/package=workflowr