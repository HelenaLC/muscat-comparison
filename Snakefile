from os.path import join
import pandas as pd
import numpy as np
import json
import re

configfile: "config.yaml"

sim_ids = open(config["sim_ids"]).read().splitlines()
agg_ids = open(config["agg_ids"]).read().splitlines()

method_ids = pd.read_csv(config["method_ids"])
method_ids = method_ids.set_index(method_ids["id"])

# put together sim_data & result directories
sim_dirs = []
res_dirs = []
fig_dirs = []
for sim_id in sim_ids:
	sim_pars = json.loads(open(join(config["sim_pars"], sim_id + ".json")).read())
	run_pars = json.loads(open(join("config", sim_id + ".json")).read())
	sim_dirs.append(\
		expand(join(config["sim_data"], "{data_id}", sim_id + ";simrep={i_sim}.rds"),\
			data_id = config["data_ids"], i_sim = range(1, sim_pars["n_reps"][0] + 1)))
	res_dirs.append(\
		expand(join(config["results"], "{data_id}", "{sim_id}",\
			"{method_id};n_genes={g}.n_cells={c}.n_samples={s};simrep={i_sim}.runrep={i_run}.rds"),\
			data_id = config["data_ids"], sim_id = sim_id, method_id = method_ids.id,\
			g = run_pars["n_genes"],\
			c = run_pars["n_cells"],\
			s = run_pars["n_samples"],\
			i_sim = range(1, sim_pars["n_reps"][0] + 1),\
			i_run = range(1, run_pars["n_reps"][0] + 1)))
	fig_dirs.append(\
		expand(join(config["figures"], "{data_id}", "fdrtpr_{sim_id};simrep={i_sim}.pdf"),\
			data_id = config["data_ids"], sim_id = sim_id,\
			i_sim = range(1, sim_pars["n_reps"][0] + 1)))
sim_dirs = sum(sim_dirs, [])
res_dirs = sum(res_dirs, [])
fig_dirs = sum(fig_dirs, [])
fig_dirs = filter(re.compile("(d[a-z][0-9]+;)|(s[0-9];)").search, fig_dirs)

rule all:
	input:	expand(join(config["raw_data"], "{data_id}-sce.rds"), data_id = config["data_ids"]),
     		expand(join(config["raw_data"], "{data_id}-sim.rds"), data_id = config["data_ids"]),
     		expand(join(config["results"],  "{data_id}_qc.html"), data_id = config["data_ids"]),
    		sim_dirs, res_dirs, fig_dirs,
    		expand(join(config["figures"], "{data_id}", "nill.pdf"), data_id = config["data_ids"]),
    		expand(join(config["figures"], "{data_id}", "upset_{sim_id}.pdf"),\
    			data_id = config["data_ids"],\
    			sim_id = ["de10", "dp10", "dm10", "db10"]),
    		expand(join(config["figures"], "{data_id}", "perf_vs_ncells_points.pdf"), data_id = config["data_ids"]),
    		expand(join(config["figures"], "{data_id}", "perf_vs_ncells_curves.pdf"), data_id = config["data_ids"]),
    		expand(join(config["figures"], "{data_id}", "perf_vs_nsamples_box.pdf"), data_id = config["data_ids"]),
    		expand(join(config["figures"], "{data_id}", "perf_split_{sim_id}.pdf"),\
    			data_id = config["data_ids"],\
    			sim_id = ["de10", "dp10", "dm10", "db10"]),
    		expand(join(config["figures"], "{data_id}", "padj_loc_vs_glb_{sim_id}.pdf"),\
    			data_id = config["data_ids"],\
    			sim_id = ["de10", "dp10", "dm10", "db10"])

# set data directories & preprocessing parameters ------------------------------
rule data_pars:
	input: 	join(config["scripts"], "data_pars.R")
	output:	expand(join(config["data_pars"], "{data_id}.rds"), data_id = config["data_ids"])
	script: "{input}"

# load data, preprocess & format into Bioconductor SCE -------------------------
rule data_sce:
	input:	script = join(config["scripts"], "data_sce.R"),
			data_pars = join(config["data_pars"], "{data_id}.rds")
	output:	join(config["raw_data"], "{data_id}-sce.rds")
	script:	"{input.script}"

# for simulation, keep only group 1 samples & fit NB ---------------------------
rule data_sim:
	input:	script = join(config["scripts"], "data_sim.R"),
			data_pars = join(config["data_pars"], "{data_id}.rds"),
			dir_data = join(config["raw_data"], "{data_id}-sce.rds")
	output:	join(config["raw_data"], "{data_id}-sim.rds")
	script:	"{input.script}"

# simulation & method parameters -----------------------------------------------
rule sim_pars:
	input:  join(config["scripts"], "sim_pars.R")
	output: config["sim_ids"]
	script: "{input}"

rule meth_pars:
    input:  script = join(config["scripts"], "meth_pars.R")
    script: "{input.script}"

# simulate data ----------------------------------------------------------------
rule sim_data:
	priority: 100
	input:  script = join(config["scripts"], "sim_data.R"),
    		data = join(config["raw_data"], "{data_id}-sim.rds"),
    		sim_pars = join(config["sim_pars"], "{sim_id}.json")
	output: sim_data = join(config["sim_data"], "{data_id}", "{sim_id};simrep={i_sim}.rds")
	script: "{input.script}"

# run DS analysis --------------------------------------------------------------
rule run_method:
	input:  script = join(config["scripts"], "run_method.R"),
			config = join("config", "{sim_id}.json"),
			sim_data = join(config["sim_data"], "{data_id}", "{sim_id};simrep={i_sim}.rds"),
			apply_fun = lambda wc: join(config["scripts"], "apply_" + method_ids.loc[wc.method_id, "type"] + ".R"),
			method_pars = lambda wc: join("metadata", method_ids.loc[wc.method_id, "type"] + "_pars.csv")
	output: res = join(config["results"], "{data_id}", "{sim_id}",\
				"{method_id};n_genes={g}.n_cells={c}.n_samples={s};simrep={i_sim}.runrep={i_run}.rds")
	script: "{input.script}"

# plot upsets
rule plot_upsets:
	wildcard_constraints: sim_id = "d[a-z][0-9]+"
	input:	script = join(config["scripts"], "plot_upsets.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.sim_id + ";").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.sim_id + "/").search, res_dirs)
	output:	fig = join(config["figures"], "{data_id}", "upset_{sim_id}.pdf")
	script:	"{input.script}"



# plot performance -------------------------------------------------------------
rule plot_perf:
	priority: 98
	input:	script = join(config["scripts"], "plot_perf.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.sim_id + ";simrep=" + wc.i_sim).search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.sim_id + "/.*simrep=" + wc.i_sim).search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "fdrtpr_{sim_id};simrep={i_sim}.pdf")
	script: "{input.script}"

# plot performance split by logFC-group ----------------------------------------
rule plot_perf_split:
	priority: 98
	input:	script = join(config["scripts"], "plot_perf_split.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.sim_id + ";").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.sim_id + "/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "perf_split_{sim_id}.pdf")
	script: "{input.script}"

# plot locally vs. globally adjusted p-values ----------------------------------
rule plot_padj_loc_vs_glb:
	input:	script = join(config["scripts"], "plot_padj_loc_vs_glb.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.sim_id + ";").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.sim_id + "/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "padj_loc_vs_glb_{sim_id}.pdf")
	script: "{input.script}"

# plot performance -------------------------------------------------------------
rule plot_perf_vs_ncells:
	wildcard_constraints: sim_id = "de10_nc"
	input:	script = join(config["scripts"], "plot_perf_vs_ncells.R"),
			utils = config["utils"],
			sim = filter(re.compile("de10_nc;").search, sim_dirs),
			res = filter(re.compile("de10_nc/").search, res_dirs)
	output: fig1 = join(config["figures"], "{data_id}", "perf_vs_ncells_points.pdf"),
			fig2 = join(config["figures"], "{data_id}", "perf_vs_ncells_curves.pdf")
	script: "{input.script}"

rule plot_perf_vs_nsamples:
	input:	script = join(config["scripts"], "plot_perf_vs_nsamples_box.R"),
			utils = config["utils"],
			sim = filter(re.compile("s[0-9]+" + ";").search, sim_dirs),
			res = filter(re.compile("s[0-9]+" + "/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "perf_vs_nsamples_box.pdf")
	script: "{input.script}"

# rule plot_time_vs_ngenes:
# 	wildcard_constraints: sim_id = "nill_ng"
# 	input:	script = join(config["scripts"], "plot_time_vs_ngenes.R"),
# 			utils = config["utils"],
# 			sim_data = join(config["sim_data"], "{data_id}", "{sim_id}.rds"),
# 			res = join(config["results"], "{data_id}", "{sim_id}")
# 	output: fig = join(config["figures"], "{data_id}", "{sim_id}.pdf")
# 	script: "{input.script}"

# null simulations -------------------------------------------------------------
rule plot_null:
	priority: 98
	input:	script = join(config["scripts"], "plot_null.R"),
			utils = config["utils"],
			res = lambda wc: filter(re.compile("nill/").search, res_dirs),
	output:	fig = join(config["figures"], "{data_id}", "nill.pdf")
	script: "{input.script}"

# # countsimQC ---------------------------------------------------------------------
# rule sim_qc:
# 	input:	script = join(config["scripts"], "countsimQC.R"),
# 			data = join(config["raw_data"], "{data_id}-sim.rds")
# 	output:	res = join(config["results"], "{data_id}_qc.html")
# 	script:	"{input.script}"

# # # runtimes ---------------------------------------------------------------------
# # rule time_methods:
# # 	input:  script = join(config["scripts"], "time_methods.R"),
# # 			dir_utils = join(config["scripts"], "utils.R"),
# # 			dir_MAST = join(config["scripts"], "MAST.R"),
# # 			dir_MM = join("scDEAmm"),
# # 			dir_data = join("data", "kang_ctrl.rds")
# # 	output: dir_res = join(config["results"], "times_methods.csv")
# # 	params: n_reps = 2, n_genes = [500, 1000, 2000, 4000, 8000],
# # 	 		methods = ["edgeR", "limma-voom", "limma-trend",\
# # 	 				   "MM-dream", "MAST", "MAST-dr"]
# # 	script: "{input.script}"

# # rule time_agg_cells:
# #     input:  script = join(config["scripts"], "time_agg_cells.R"),
# #     		dir_utils = join(config["scripts"], "utils.R"),
# #     		dir_data = join(config["data"], "kang_ctrl.rds"),
# #     		dir_agg_pars = config["dir_agg_pars"]
# #     output: dir_res = join(config["results"], "times_agg_cells.csv")
# #     params: n_reps = 10, n_cells = [50, 100, 200, 400, 800]
# #     script: "{input.script}"

# # rule time_agg_genes:
# #     input:  script = join(config["scripts"], "time_agg_genes.R"),
# #     		dir_utils = join(config["scripts"], "utils.R"),
# #     		dir_data = join("data", "kang_ctrl.rds"),
# #     		dir_agg_pars = config["dir_agg_pars"]
# #     output: dir_res = join(config["results"], "times_agg_genes.csv")
# #     params: n_reps = 10, n_genes = [500, 1000, 2000, 4000, 8000]
# #     script: "{input.script}"

# # rule plot_times_agg_cells:
# # 	input:  script = join(config["scripts"], "plot_times_agg_cells.R"),
# # 			dir_utils = join(config["scripts"], "utils.R"),
# # 			dir_res = join(config["results"], "times_agg_cells.csv")
# # 	output:	dir_fig = join(config["figures"], "times_agg_cells.pdf")
# # 	script: "{input.script}"

# # rule plot_times_agg_genes:
# # 	input:  script = join(config["scripts"], "plot_times_agg_genes.R"),
# # 			dir_utils = join(config["scripts"], "utils.R"),
# # 			dir_res = join(config["results"], "times_agg_genes.csv")
# # 	output:	dir_fig = join(config["figures"], "times_agg_genes.pdf")
# # 	script: "{input.script}"

# # rule plot_times_methods:
# # 	input:  script = join(config["scripts"], "plot_times_methods.R"),
# # 			dir_utils = join(config["scripts"], "utils.R"),
# # 			dir_res = join(config["results"], "times_methods.csv")
# # 	output:	dir_fig = join(config["figures"], "times_methods.pdf")
# # 	script: "{input.script}"
