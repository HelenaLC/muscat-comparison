onstart:
	shell("Rscript scripts/data_pars.R")
	shell("Rscript scripts/sim_pars.R")
	shell("Rscript scripts/meth_pars.R")

from os.path import join
import pandas as pd
import numpy as np
import json
import re

configfile: "config.yaml"

sim_ids = open(config["sim_ids"]).read().splitlines()
method_ids = pd.read_csv(config["method_ids"])
method_ids = method_ids.set_index(method_ids["id"])

# put together sim_data & result directories
sim_dirs = []
res_dirs = []
fig_dirs = []
prf_dirs = []
for dat_id in config["dat_ids"]:
	for sim_id in sim_ids:
		sim_pars = json.loads(open(join(config["sim_pars"], sim_id + ".json")).read())
		run_pars = json.loads(open(join("config", dat_id + ";" + sim_id + ".json")).read())
		if (run_pars == ['norun']): continue
		sim_is = range(1, sim_pars["n_reps"][0] + 1)
		run_is = range(1, run_pars["n_reps"][0] + 1)
		sim_dirs.append(\
			expand(join(config["sim_data"], "{data_id}", sim_id + ";simrep={i_sim}.rds"),\
				data_id = dat_id, i_sim = sim_is))
		res_dirs.append(\
			expand(join(config["results"], "{data_id}", "{sim_id}",\
				"{method_id};n_genes={g}.n_cells={c}.n_samples={s};simrep={i_sim}.runrep={i_run}.rds"),\
				data_id = dat_id, sim_id = sim_id, method_id = method_ids.id,\
				g = run_pars["n_genes"],\
				c = run_pars["n_cells"],\
				s = run_pars["n_samples"],\
				i_sim = sim_is,\
				i_run = run_is))
		prf_dirs.append(expand(
			join(config["results"], "{data_id}", "cobra", "{sim_id};{method_id};si={i_sim};ri={i_run}.rds"),\
			data_id = dat_id, sim_id = sim_id, method_id = method_ids.id, i_sim = sim_is, i_run = run_is))
		fig_dirs.append(\
			expand(join(config["figures"], "{data_id}", "{sim_id}", "fdrtpr_{sim_id};simrep={i_sim}.pdf"),\
				data_id = dat_id, sim_id = sim_id, i_sim = sim_is))

sim_dirs = sum(sim_dirs, [])
res_dirs = sum(res_dirs, [])
prf_dirs = sum(prf_dirs, [])
fig_dirs = sum(fig_dirs, [])
fig_dirs = filter(re.compile("(d[a-z][0-9]+;)|(mixed_ss;)").search, fig_dirs)

rule all:
	input:	expand(join(config["raw_data"], "{data_id}-sce.rds"), data_id = config["dat_ids"]),
     		expand(join(config["raw_data"], "{data_id}-sim.rds"), data_id = config["dat_ids"]),
     		#expand(join(config["results"],  "{data_id}_qc.html"), data_id = config["dat_ids"]),
     		#expand(join(config["figures"],  "{data_id}", "pb_mean_disp.pdf"), data_id = config["dat_ids"]),
     		#expand(join(config["figures"],  "{data_id}", "pbMDS.pdf"), data_id = config["dat_ids"]),
    		sim_dirs, res_dirs, prf_dirs, #fig_dirs,
    		expand(join(config["figures"], "{data_id}", "nill.pdf"), data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "upsets_split_by_cat.pdf"),\
    		#	data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "perf_vs_ncells_points.pdf"), data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "perf_vs_ncells_curves.pdf"), data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "perf_vs_nsamples_box.pdf"), data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "padj_local_vs_global.pdf"), data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "perf_split_by_cat.pdf"),    data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "perf_split_by_lfc.pdf"),    data_id = config["dat_ids"]),
    		#expand(join(config["figures"], "{data_id}", "sim_vs_est_logFC.pdf"),     data_id = config["dat_ids"]),
    		#expand(join(config["results"], "{data_id}", "{sim_id};concordance.rds"),\
    		#	data_id = config["dat_ids"], sim_id = ["de10", "dp10", "dm10"]),
    		#expand(join(config["figures"], "{data_id}", "tree_{sim_id}.pdf"),\
    		#	data_id = config["dat_ids"], sim_id = ["de10", "dp10", "dm10"])

# countsimQC -------------------------------------------------------------------
# rule sim_qc:
# 	input:	script = join(config["scripts"], "countsimQC.R"),
# 			data = join(config["raw_data"], "{data_id}-sim.rds")
# 	output:	res = join(config["results"], "{data_id}_qc.html")
# 	script:	"{input.script}"

# load data, preprocess & format into Bioconductor SCE -------------------------
rule data_sce:
	priority: 100
	input:	script = join(config["scripts"], "data_sce.R"),
			data_pars = join(config["data_pars"], "{data_id}.rds")
	output:	join(config["raw_data"], "{data_id}-sce.rds")
	script:	"{input.script}"

# for simulation, keep only group 1 samples & fit NB ---------------------------
rule data_sim:
	priority: 100
	input:	script = join(config["scripts"], "data_sim.R"),
			data_pars = join(config["data_pars"], "{data_id}.rds"),
			dir_data = join(config["raw_data"], "{data_id}-sce.rds")
	output:	join(config["raw_data"], "{data_id}-sim.rds")
	script:	"{input.script}"

# pb-level mean-disperion ------------------------------------------------------
rule pb_mean_disp:
	priority: 99
	input:	script = join(config["plotting"], "plot_pb_mean_disp.R"),
			utils = config["utils"],
			data = join(config["raw_data"], "{data_id}-sim.rds")
	output:	gg = join(config["figures"], "{data_id}", "pb_mean_disp.rds"),
			fig = join(config["figures"], "{data_id}", "pb_mean_disp.pdf")
	script:	"{input.script}"

# # pb-level MDS plots -----------------------------------------------------------
# rule plot_pbMDS:
# 	priority: 99
# 	input:	script = join(config["scripts"], "plot_pbMDS.R"),
# 			utils = config["utils"],
# 			data = join(config["raw_data"], "{data_id}-sim.rds")
# 	output:	fig = join(config["figures"], "{data_id}", "pbMDS.pdf")
# 	script:	"{input.script}"

# data simulation --------------------------------------------------------------
rule sim_data:
	priority: 99
	input:  script = join(config["scripts"], "sim_data.R"),
    		data = join(config["raw_data"], "{data_id}-sim.rds"),
    		sim_pars = join(config["sim_pars"], "{sim_id}.json")
	output: sim_data = join(config["sim_data"], "{data_id}", "{sim_id};simrep={i_sim}.rds")
	script: "{input.script}"

# DS analysis ------------------------------------------------------------------
rule run_method:
	input:  script = join(config["scripts"], "run_method.R"),
			config = join("config", "{data_id};{sim_id}.json"),
			sim_data = join(config["sim_data"], "{data_id}", "{sim_id};simrep={i_sim}.rds"),
			apply_fun = lambda wc: join(config["scripts"], "apply_" + method_ids.loc[wc.method_id, "type"] + ".R"),
			method_pars = lambda wc: join(config["meth_pars"], method_ids.loc[wc.method_id, "type"] + "_pars.csv")
	output: res = join(config["results"], "{data_id}", "{sim_id}",\
				"{method_id};n_genes={g}.n_cells={c}.n_samples={s};simrep={i_sim}.runrep={i_run}.rds")
	script: "{input.script}"

rule calc_perf:
	input:  script = join(config["scripts"], "calc_perf.R"),
			sce = join(config["sim_data"], "{data_id}", "{sim_id};simrep={i_sim}.rds"),
			res = join(config["results"], "{data_id}", "{sim_id}", "{method_id};n_genes={g}.n_cells={c}.n_samples={s};simrep={i_sim}.runrep={i_run}.rds")
	output: prf = config["results"], "{data_id}", "cobra", "{sim_id};{method_id};si={i_sim};ri={i_run}.rds")
	script: "{input.script}"

# plot upsets
rule plot_upsets:
	wildcard_constraints: sim_id = "d[a-z][0-9]+"
	input:	script = join(config["plotting"], "plot_upsets.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile("d[a-z][0-9]+;").search, sim_dirs),
			res = lambda wc: filter(re.compile("d[a-z][0-9]+/").search, res_dirs)
	output:	fig = join(config["figures"], "{data_id}", "upsets_split_by_cat.pdf")
	script:	"{input.script}"

# plot performance -------------------------------------------------------------
rule plot_perf:
	input:	script = join(config["plotting"], "plot_perf.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.sim_id + ";simrep=" + wc.i_sim).search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/" + wc.sim_id + "/.*simrep=" + wc.i_sim).search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "{sim_id}", "fdrtpr_{sim_id};simrep={i_sim}.pdf")
	script: "{input.script}"

# plot locally vs. globally adjusted p-values ----------------------------------
rule plot_padj_loc_vs_glb:
	#wildcard_constraints: sim_id = "d[a-z][0-9]+$"
	input:	script = join(config["plotting"], "plot_padj_local_vs_global.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+;").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "padj_local_vs_global.pdf"),
			gg = join(config["figures"], "{data_id}", "padj_local_vs_global.rds")
	script: "{input.script}"

# plot performance split by gene catergoy --------------------------------------
rule plot_perf_split_by_cat:
	input:	script = join(config["plotting"], "plot_perf_split_by_cat.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+;").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "perf_split_by_cat.pdf"),
			gg = join(config["figures"], "{data_id}", "perf_split_by_cat.rds")
	script: "{input.script}"

# plot performance split by logFC-group ----------------------------------------
rule plot_perf_split_by_lfc:
	input:	script = join(config["plotting"], "plot_perf_split_by_lfc.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+;").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "perf_split_by_lfc.pdf")
			#gg = join(config["figures"], "{data_id}", "perf_split_by_cat.rds")
	script: "{input.script}"

# plot performance -------------------------------------------------------------
rule plot_perf_vs_ncells:
	input:	script = join(config["plotting"], "plot_perf_vs_ncells.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.data_id + "/de10_nc;").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/de10_nc/").search, res_dirs)
	output: fig1 = join(config["figures"], "{data_id}", "perf_vs_ncells_points.pdf"),
			fig2 = join(config["figures"], "{data_id}", "perf_vs_ncells_curves.pdf"),
			gg = join(config["figures"], "{data_id}", "perf_vs_ncells.rds")
	script: "{input.script}"

rule plot_perf_vs_nsamples:
	input:	script = join(config["plotting"], "plot_perf_vs_nsamples.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.data_id + "/de10_ns;").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/de10_ns/").search, res_dirs)
	output: fig1 = join(config["figures"], "{data_id}", "perf_vs_nsamples_points.pdf"),
			fig2 = join(config["figures"], "{data_id}", "perf_vs_nsamples_curves.pdf"),
			gg = join(config["figures"], "{data_id}", "perf_vs_nsamples.rds")
	script: "{input.script}"

rule plot_sim_vs_est_logFC:
	input:	script = join(config["plotting"], "plot_sim_vs_est_logFC.R"),
			utils = config["utils"],
			sim = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+;").search, sim_dirs),
			res = lambda wc: filter(re.compile(wc.data_id + "/d[a-z][0-9]+/").search, res_dirs)
	output: fig = join(config["figures"], "{data_id}", "sim_vs_est_logFC.pdf")
	script: "{input.script}"

# # rule plot_time_vs_ngenes:
# # 	wildcard_constraints: sim_id = "nill_ng"
# # 	input:	script = join(config["plotting"], "plot_time_vs_ngenes.R"),
# # 			utils = config["utils"],
# # 			sim_data = join(config["sim_data"], "{data_id}", "{sim_id}.rds"),
# # 			res = join(config["results"], "{data_id}", "{sim_id}")
# # 	output: fig = join(config["figures"], "{data_id}", "{sim_id}.pdf")
# # 	script: "{input.script}"

# null simulations -------------------------------------------------------------
rule plot_null:
	priority: 98
	input:	script = join(config["plotting"], "plot_null.R"),
			utils = config["utils"],
			res = lambda wc: filter(re.compile(wc.data_id + "/nill/").search, res_dirs),
	output:	fig = join(config["figures"], "{data_id}", "nill.pdf"),
			res = join(config["results"], "{data_id}", "nill_df.rds")
	script: "{input.script}"

# cross-method concordance ------------------------------------------------------
# rule calc_conc:
# 	priority: 98
# 	input:	script = join(config["scripts"], "calc_conc.R"),
# 			utils = config["utils"],
# 			sim = lambda wc: filter(re.compile(wc.sim_id + ";").search, sim_dirs),
# 			res = lambda wc: filter(re.compile(wc.sim_id + "/").search, res_dirs)
# 	output:	res = join(config["results"], "{data_id}", "{sim_id};concordance.rds")
# 	script: "{input.script}"

# rule plot_conc:
# 	priority: 98
# 	input:	script = join(config["plotting"], "plot_conc.R"),
# 			utils = config["utils"],
# 			res = join(config["results"], "{data_id}", "{sim_id};concordance.rds")
# 	output:	fig = join(config["figures"], "{data_id}", "tree_{sim_id}.pdf")
# 	script: "{input.script}"
