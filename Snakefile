# did: reference data ID
# sid: simulation ID
# mid: method ID

# i: simulation replicate
# j: run replicate

# g: # genes
# c: # cells
# k: # clusters
# s: # samples

onstart:
	shell("Rscript scripts/sim_pars.R")
	shell("Rscript scripts/run_pars.R")
	shell("Rscript scripts/meth_pars.R")

configfile: "config.yaml"

from os.path import join
import pandas as pd
import json
import re

sids = json.loads(open(config["sids"]).read())
mids = pd.read_csv(config["mids"])
mids = mids.set_index(mids["id"])

sim_dirs = []
res_dirs = []
fig_dirs = []
for sid in sids:
	sim_pars = json.loads(open(join(config["sim_pars"], sid + ".json")).read())
	for did in config["dids"]:
		run_pars = json.loads(open(join(config["run_pars"], did + ";" + sid + ".json")).read())
		if run_pars is None: continue
		sim_dirs.append(expand(join(\
			config["sim_data"], "{did}", "{sid};{i}.rds"),\
			did = did, sid = sid, i = range(1, sim_pars["nr"][0] + 1)))
		fig_dirs.append(expand(join(\
			config["figures"], "{did}", "tprfdr", "{sid};{i}.pdf"),\
			did = did, sid = sid, i = range(1, sim_pars["nr"][0] + 1)))
		res_dirs.append(expand(join(\
			config["results"], "{did}", "{sid};{i};{mid};{j};g{g};c{c};k{k};s{s}.rds"),
			did = did, sid = sid, mid = mids.id,\
			i = range(1, sim_pars["nr"][0] + 1),\
			j = range(1, run_pars["nr"][0] + 1),\
			g = run_pars["ng"], c = run_pars["nc"],\
			k = run_pars["nk"], s = run_pars["ns"]))

sim_dirs = sum(sim_dirs, [])
res_dirs = sum(res_dirs, [])
fig_dirs = sum(fig_dirs, [])
fig_dirs = filter(re.compile("d[a-z][0-9]+;").search, fig_dirs)

rule all: 
	input:	sim_dirs, res_dirs, fig_dirs,
			#expand(join(config["figures"], "{did}_qc.html"), did = config["dids"]),
			expand(join(config["figures"], "{did}", "{nms}.pdf"), did = config["dids"],\
				nms = ["pb_mean_disp", "perf_by_cat", "perf_by_ss"]),
			expand(join(config["figures"], "{did}", "upset.pdf"), did = config["dids"]),
			expand(join(config["figures"], "{did}", "perf_by_cat_{padj}.{ext}"),\
				did = config["dids"], padj = ["loc", "glb"], ext = ["rds", "pdf"]),
			expand(join(config["figures"], "{did}", "{nms}.{ext}"), did = config["dids"],\
				nms = ["nill", "sim_vs_est_lfc"], ext = ["rds", "pdf"]),
			expand(join(config["figures"], "{did}", "perf_by_n{x}.{ext}"),\
				did = config["dids"], x = "c", ext = ["rds", "pdf"]),
			expand(join(config["figures"], "{did}", "perf_by_n{x}.{ext}"),\
				did = "kang", x = "s", ext = ["rds", "pdf"]),
			expand(join(config["figures"], "{did}", "perf_by_expr_{padj}.{ext}"),\
				did = config["dids"], padj = ["loc", "glb"], ext = ["rds", "pdf"])

rule sim_qc:
	input:	script = join(config["scripts"], "sim_qc.R"),
			sce = lambda wc: join(config["raw_data"], wc.did + ".rds")
	output: html = join(config["figures"], "{did}_qc.html")
	script:	"{input.script}"

rule sim_data:
	priority: 100
	input:  script = join(config["scripts"], "sim_data.R"),
			sce = lambda wc: join(config["raw_data"], wc.did + ".rds"),
			sim_pars = join(config["sim_pars"], "{sid}.json")
	output: sim = join(config["sim_data"], "{did}", "{sid};{i}.rds")
	script: "{input.script}"

rule run_meth:
	priority: 99
	input:	script = join(config["scripts"], "run_meth.R"),
			sim = join(config["sim_data"], "{did}", "{sid};{i}.rds"),
			run_pars = join(config["run_pars"], "{did};{sid}.json"),
			meth_pars = join(config["meth_pars"], "{mid}.json"),
			fun = lambda wc: join(config["scripts"], "apply_" + mids.loc[wc.mid, "type"] + ".R")
	output: res = join(config["results"], "{did}", "{sid};{i};{mid};{j};g{g};c{c};k{k};s{s}.rds")
	script:	"{input.script}"

rule plot_pb_mean_disp:
	input:	config["utils"],
			script = join(config["scripts"], "plot_pb_mean_disp.R"),
			sce = lambda wc: join(config["raw_data"], wc.did + ".rds"),
	output: ggp = join(config["figures"], "{did}", "pb_mean_disp.rds"),
			fig = join(config["figures"], "{did}", "pb_mean_disp.pdf")
	script: "{input.script}"

rule plot_null:
	input:	config["utils"],
			script = join(config["scripts"], "plot_null.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/nill;.*").search, res_dirs)
	output: fig = join(config["figures"], "{did}", "null.pdf")
	script: "{input.script}"

rule plot_tprfdr:
	input:	config["utils"],
			script = join(config["scripts"], "plot_tprfdr.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did + "/" +\
				wc.sid + ";" + wc.i + ";.*").search, res_dirs)
	output: fig = join(config["figures"], "{did}", "tprfdr", "{sid};{i}.pdf")
	script: "{input.script}"

rule plot_perf_by_cat:
	input:	config["utils"],
			script = join(config["scripts"], "plot_perf_by_cat.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/d[a-z][0-9]+;.*").search, res_dirs)
	output: ggp = join(config["figures"], "{did}", "perf_by_cat_{padj}.rds"),
			fig = join(config["figures"], "{did}", "perf_by_cat_{padj}.pdf")
	script: "{input.script}"

rule plot_perf_by_nx:
	input:	config["utils"],
			script = join(config["scripts"], "plot_perf_by_nx.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/ds10_n" + wc.x + ";.*").search, res_dirs)
	output: ggp = join(config["figures"], "{did}", "perf_by_n{x}.rds"),
			fig = join(config["figures"], "{did}", "perf_by_n{x}.pdf")
	script: "{input.script}"

rule plot_perf_by_ss:
	input:	config["utils"],
			script = join(config["scripts"], "plot_perf_by_ss.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/ds10_ss").search, res_dirs)
	output: fig = join(config["figures"], "{did}", "perf_by_ss.pdf")
	script: "{input.script}"

rule plot_perf_by_expr:
	input:	config["utils"],
			script = join(config["scripts"], "plot_perf_by_expr.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/ds10;.*").search, res_dirs)
	output: ggp = join(config["figures"], "{did}", "perf_by_expr_{padj}.rds"),
			fig = join(config["figures"], "{did}", "perf_by_expr_{padj}.pdf")
	script: "{input.script}"

rule plot_upset:
	input:	config["utils"],
			script = join(config["scripts"], "plot_upset.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/d[a-z][0-9]+;.*").search, res_dirs)
	output: fig = join(config["figures"], "{did}", "upset.pdf")
	script: "{input.script}"

rule plot_lfc:
	input:	config["utils"],
			script = join(config["scripts"], "plot_lfc.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				"/d[a-z]10;.*").search, res_dirs)
	output:	ggp = join(config["figures"], "{did}", "sim_vs_est_lfc.rds"),
			fig = join(config["figures"], "{did}", "sim_vs_est_lfc.pdf")
	script: "{input.script}"

rule fig_perf_by_cat:
	input:	config["utils"],
			script = join(config["scripts"], "fig_perf_by_cat.R"),
			gg1 = lambda wc: join(config["figures"] + "/" + wc.did + "/perf_by_cat_loc.rds"),
			gg2 = lambda wc: join(config["figures"] + "/" + wc.did + "/perf_by_cat_glb.rds")
	output: fig = join(config["figures"], "{did}", "perf_by_cat.pdf")
	script: "{input.script}"







