# did: reference dataset ID
# sid: simulation ID
# mid: method ID

# i: simulation replicate
# j: method run replicate

# g: # genes
# c: # cells
# k: # clusters
# s: # samples

onstart:
	shell("Rscript scripts/sim_pars.R")  # simulation parameters
	shell("Rscript scripts/run_pars.R")  # run mode parameters
	shell("Rscript scripts/meth_pars.R") # method parameters

configfile: "config.yaml"

from os.path import join
import pandas as pd
import json
import re

os.environ['KMP_DUPLICATE_LIB_OK']='True'

sids = json.loads(open(config["sids"]).read())
mids = pd.read_csv(config["mids"])
mids = mids.set_index(mids["id"])

sim_dirs = []
res_dirs = []
for sid in sids:
	sim_pars = json.loads(open(join(config["sim_pars"], sid + ".json")).read())
	for did in config["dids"]:
		run_pars = json.loads(open(join(config["run_pars"], did + ";" + sid + ".json")).read())
		if run_pars is None: continue
		sim_dirs.append(expand(join(\
			config["sim_data"], "{did};{sid};{i}.rds"),\
			did = did, sid = sid, i = range(1, sim_pars["nr"][0] + 1)))
		res_dirs.append(expand(join(\
			config["results"], "{did};{sid};{i};{mid};{j};g{g};c{c};k{k};s{s}.rds"),
			did = did, sid = sid, mid = mids.id,\
			i = range(1, sim_pars["nr"][0] + 1),\
			j = range(1, run_pars["nr"][0] + 1),\
			g = run_pars["ng"], c = run_pars["nc"],\
			k = run_pars["nk"], s = run_pars["ns"]))

sim_dirs = sum(sim_dirs, [])
res_dirs = sum(res_dirs, [])

rule all: 
	input:	sim_dirs, res_dirs,
			#expand(join(config["raw_data"], "sce_{did}.rds"), did = config["dids"]),
			#expand(join(config["raw_data"], "ref_{did}.rds"), did = config["dids"]),
		# null simulation p-value distributions 
			expand(join(config["figures"], "{did}-null.{ext}"),\
				did = config["dids"], ext = ["rds", "pdf"]),
		#	expand(join(config["figures"], "{did}_qc.html"), did = config["dids"]),
			expand(join(config["figures"], "{did}-pb_mean_disp.pdf"), did = config["dids"]),
			expand(join(config["figures"], "{did}-upset.pdf"), did = config["dids"]),
		# TPR-FDR stratified by DD category using locally/globally adjusted p-values
			expand(join(config["figures"], "{did}-perf_by_cat_{padj}.{ext}"),\
				did = config["dids"], padj = ["loc", "glb"], ext = ["rds", "pdf"])
		#	expand(join(config["figures"], "{did}", "{nms}.{ext}"), did = config["dids"],\
		#	nms = ["null", "perf_by_ss", "sim_vs_est_lfc"], ext = ["rds", "pdf"]),
			#expand(join(config["figures"], "{did}", "perf_by_n{x}.{ext}"),\
		#		did = config["dids"], x = "c", ext = ["rds", "pdf"]),
			#expand(join(config["figures"], "{did}", "perf_by_n{x}.{ext}"),\
		#		did = "kang", x = "s", ext = ["rds", "pdf"]),
		#	expand(join(config["figures"], "{did}", "perf_by_expr_{padj}.{ext}"),\
		#		did = config["dids"], padj = ["loc", "glb"], ext = ["rds", "pdf"]),
		# method runtimes versus no. cells/genes
		#	expand(join(config["figures"], "{did}", "runtimes.pdf"), did = ["kang"]),
		#	#expand(join(config["results"], "lps", "{mid}.rds"), mid = mids.id)

rule prep_sce:
	input:  script = join(config["scripts"], "prep_{did}.R"),
    		sce = join(config["raw_data"], "sce0_{did}.rds")
	output:	sce = join(config["raw_data"], "sce_{did}.rds")
	script:	"scripts/prep_kang.R"
	# shell:	'''R CMD BATCH --no-restore --no-save "--args\
	# 	input_sce='{input.input_sce}' output_sce='{output.output_sce}'"\
	# 	{input.script} {log}'''

rule prep_sim:
	input:	script = join(config["scripts"], "prep_sim.R"),
			sce = join(config["raw_data"], "sce_{did}.rds")
	output:	sce = join(config["raw_data"], "ref_{did}.rds")
	script:	"scripts/prep_sim.R"
	# shell:	'''R CMD BATCH --no-restore --no-save "--args\
	# 	input_sce='{input.input_sce}' output_sce='{output.output_sce}'"\
	# 	{input.script} {log}'''

# rule sim_qc:
# 	priority: 98
# 	input:	script = join(config["scripts"], "sim_qc.R"),
# 			sce = lambda wc: join(config["raw_data"], "ref_" + wc.did + ".rds")
# 	output: html = join(config["figures"], "{did}_qc.html")
# 	script:	"{input.script}"

rule sim_data:
	input:  script = join(config["scripts"], "sim_data.R"),
			sim_pars = join(config["sim_pars"], "{sid}.json"),
			sce = join(config["raw_data"], "ref_{did}.rds")
	output: sce = join(config["sim_data"], "{did};{sid};{i}.rds")
	script:	"scripts/sim_data.R"

rule run_meth:
	input:	script = join(config["scripts"], "run_meth.R"),
			sim = join(config["sim_data"], "{did};{sid};{i}.rds"),
			meth_pars = join(config["meth_pars"], "{mid}.json"),
			run_pars = join(config["run_pars"], "{did};{sid}.json"),
			fun = lambda wc: join(config["scripts"], "apply_" + mids.loc[wc.mid, "type"] + ".R")
	output: res = join(config["results"], "{did};{sid};{i};{mid};{j};g{g};c{c};k{k};s{s}.rds")
	script:	"scripts/run_meth.R"

# pseudobulk-level mean-dispersion scatters
rule plot_pb_mean_disp:
	input:	config["utils"],
			script = join(config["scripts"], "plot_pb_mean_disp.R"),
			sce = join(config["raw_data"], "ref_{did}.rds"),
	output: ggp = join(config["figures"], "{did}-pb_mean_disp.rds"),
			fig = join(config["figures"], "{did}-pb_mean_disp.pdf")
	script: "scripts/plot_pb_mean_disp.R"

# null simulation p-value distributions
rule plot_null:
	input:	config["utils"],
			script = join(config["scripts"], "plot_null.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				";nill;.*").search, res_dirs)
	output: ggp = join(config["figures"], "{did}-null.rds"), 
			fig = join(config["figures"], "{did}-null.pdf")
	script:	"scripts/plot_null.R"

# TPR-FDR stratified by DD category using locally/globally adjusted p-values
rule plot_perf_by_cat:
	input:	config["utils"],
			script = join(config["scripts"], "plot_perf_by_cat.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				";d[a-z][0-9]+;.*").search, res_dirs)
	output: ggp = join(config["figures"], "{did}-perf_by_cat_{padj}.rds"),
			fig = join(config["figures"], "{did}-perf_by_cat_{padj}.pdf")
	script: "scripts/plot_perf_by_cat.R"

# TPR-FDR points stratified by #(cells)/#(samples)
# rule plot_perf_by_nx:
# 	input:	config["utils"],
# 			script = join(config["scripts"], "plot_perf_by_nx.R"),
# 			res = lambda wc: filter(re.compile(\
# 				config["results"] + "/" + wc.did +\
# 				"/de10_n" + wc.x + ";.*").search, res_dirs)
# 	output: ggp = join(config["figures"], "{did}", "perf_by_n{x}.rds"),
# 			fig = join(config["figures"], "{did}", "perf_by_n{x}.pdf")
# 	script: "{input.script}"

# TPR-FDR points stratified by #(samples)
# rule plot_perf_by_ss:
# 	input:	config["utils"],
# 			script = join(config["scripts"], "plot_perf_by_ss.R"),
# 			res = lambda wc: filter(re.compile(\
# 				config["results"] + "/" + wc.did +\
# 				"/de10_ss").search, res_dirs)
# 	output: ggp = join(config["figures"], "{did}", "perf_by_ss.rds"),
# 			fig = join(config["figures"], "{did}", "perf_by_ss.pdf")
# 	script: "{input.script}"

# TPR-FDR points stratified by expression
# rule plot_perf_by_expr:
# 	input:	config["utils"],
# 			script = join(config["scripts"], "plot_perf_by_expr.R"),
# 			res = lambda wc: filter(re.compile(\
# 				config["results"] + "/" + wc.did +\
# 				"/de10;.*").search, res_dirs)
# 	output: ggp = join(config["figures"], "{did}", "perf_by_expr_{padj}.rds"),
# 			fig = join(config["figures"], "{did}", "perf_by_expr_{padj}.pdf")
# 	script: "{input.script}"

rule plot_upset:
	input:	config["utils"],
			script = join(config["scripts"], "plot_upset.R"),
			res = lambda wc: filter(re.compile(\
				config["results"] + "/" + wc.did +\
				";d[a-z][0-9]+;.*").search, res_dirs)
	output: fig = join(config["figures"], "{did}-upset.pdf")
	script: "scripts/plot_upset.R"

# scatter plots of estimated vs. simulated logFCs
# rule plot_lfc:
# 	input:	config["utils"],
# 			script = join(config["scripts"], "plot_lfc.R"),
# 			res = lambda wc: filter(re.compile(\
# 				config["results"] + "/" + wc.did +\
# 				"/d[a-z]10;.*").search, res_dirs)
# 	output:	ggp = join(config["figures"], "{did}", "sim_vs_est_lfc.rds"),
# 			fig = join(config["figures"], "{did}", "sim_vs_est_lfc.pdf")
# 	script: "{input.script}"

# runtimes vs. #(cells)/#(genes)
# rule plot_runtimes:
# 	input: 	config["utils"],
# 			script = join(config["scripts"], "plot_runtimes.R"),
# 			res = lambda wc: filter(re.compile(\
# 				config["results"] + "/" + wc.did +\
# 				"/de10_n[g|c];*").search, res_dirs)
# 	output:	fig = join(config["figures"], "{did}", "runtimes.pdf")
# 	script:	"{input.script}"

#rule lps_run_meth:
#	input:	script = join(config["scripts"], "lps_run_meth.R"),
#			sce = join(config["raw_data"], "sce0_magl.rds"),
#			meth_pars = join(config["meth_pars"], "{mid}.json"),
#			fun = lambda wc: join(config["scripts"], "apply_" + mids.loc[wc.mid, "type"] + ".R")
#	output:	res = join(config["results"], "lps", "{mid}.rds")
#	script:	"{input.script}"

# write session info to .txt file
# rule session_info:
# 	input:	join(config["scripts"], "session_info.R")
# 	output:	txt = "session_info.txt"
# 	script:	"{input}"







