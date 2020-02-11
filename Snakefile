# did: reference dataset ID
# sid: simulation ID
# mid: method ID

# i: simulation replicate
# j: method run replicate

# g: # genes
# c: # cells
# k: # clusters
# s: # samples

import pandas as pd
import json
import re

onstart:
	shell("Rscript scripts/sim_pars.R")  # simulation parameters
	shell("Rscript scripts/run_pars.R")  # run mode parameters
	shell("Rscript scripts/meth_pars.R") # method parameters

configfile: "config.yaml"

R = config["R"]

#os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

sids = json.loads(open(config["sids"]).read())
mids = pd.read_csv(config["mids"])
mids = mids.set_index(mids["id"])

sim_dirs = []
res_dirs = []
for sid in sids:
	sim_pars = json.loads(open(config["sim_pars"] + sid + ".json").read())
	for did in config["dids"]:
		run_pars = json.loads(open(config["run_pars"] + did + "," + sid + ".json").read())
		if run_pars is None: continue
		sim_dirs.append(expand(\
			config["sim_data"] + "{did},{sid},{i}.rds",\
			did = did, sid = sid, i = range(1, sim_pars["nr"][0] + 1)))
		res_dirs.append(expand(\
			config["results"] + "{did},{sid},{i},{mid},{j},g{g},c{c},k{k},s{s}.rds",
			did = did, sid = sid, mid = mids.id,\
			i = range(1, sim_pars["nr"][0] + 1),\
			j = range(1, run_pars["nr"][0] + 1),\
			g = run_pars["ng"], c = run_pars["nc"],\
			k = run_pars["nk"], s = run_pars["ns"]))

sim_dirs = sum(sim_dirs, [])
res_dirs = sum(res_dirs, [])

rule all: 
	input:	sim_dirs, res_dirs,
			#expand(config["raw_data"], "sce_{did}.rds"), did = config["dids"]),
			#expand(config["raw_data"], "ref_{did}.rds"), did = config["dids"]),
		# 'countsimQC' reports
			expand(config["plots"] + "{did}-sim_qc.html", did = config["dids"]),
		# null simulation p-value distributions 
			expand(config["plots"] + "{did}-null.{ext}",\
				did = config["dids"], ext = ["rds", "pdf"]),
		# expand(config["plots"], "{did}_qc.html"), did = config["dids"]),
			expand(config["plots"] + "{did}-pb_mean_disp.pdf", did = config["dids"]),
			expand(config["plots"] + "{did}-upset.pdf", did = config["dids"]),
		# TPR-FDR stratified by DD category using locally/globally adjusted p-values
			expand(config["plots"] + "{did}-perf_by_cat_{padj}.{ext}",\
				did = config["dids"], padj = ["loc", "glb"], ext = ["rds", "pdf"]),
		# scatters of simulated vs. estimated logFCs
			expand(config["plots"] + "{did}-sim_vs_est_lfc.{ext}",\
				did = config["dids"], ext = ["rds", "pdf"]),
		# TPR-FDR stratified by no. of cells per cluster-sample
			expand(config["plots"] + "{did}-perf_by_n{x}.{ext}",\
				did = config["dids"], x = "c", ext = ["rds", "pdf"]),
		# TPR-FDR stratified by no. of replicates per group
			expand(config["plots"] + "{did}-perf_by_n{x}.{ext}",\
				did = "kang", x = "s", ext = ["rds", "pdf"]),
		# TPR-FDR stratified by magnitude of sample-size unbalancing
			expand(config["plots"] + "{did}-perf_by_{y}s.{ext}",\
				did = config["dids"], y = ["s", "g"], ext = ["rds", "pdf"]),
		# TPR-FDR stratified by expression level 
			expand(config["plots"] + "{did}-perf_by_es_{padj}.{ext}",\
				did = config["dids"], padj = ["loc", "glb"], ext = ["rds", "pdf"]),
		# method runtimes versus no. cells/genes
			expand(config["plots"] + "{did}-runtimes.pdf", did = ["kang"])
			#expand(config["results"], "lps", "{mid}.rds"), mid = mids.id)

rule prep_sce:
	priority: 100
	input:  script = config["scripts"] + "prep_{did}.R",
    		sce = config["raw_data"] + "sce0_{did}.rds"
	output:	sce = config["raw_data"] + "sce_{did}.rds"
	log:	config["logs"] + "prep_sce-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args input_sce={input.sce} output_sce={output.sce}"\
		{input.script} {log}'''

rule prep_sim:
	priority: 99
	input:	script = config["scripts"] + "prep_sim.R",
			sce = config["raw_data"] + "sce_{did}.rds"
	output:	sce = config["raw_data"] + "ref_{did}.rds"
	log:	config["logs"] + "prep_sim-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args input_sce={input.sce} output_sce={output.sce}"\
		{input.script} {log}'''

rule sim_qc:
	priority: 98
	input:	script = config["scripts"] + "sim_qc.R",
			sce = config["raw_data"] + "ref_{did}.rds"
	output: config["plots"] + "{did}-sim_qc.html"
	log:	config["logs"] + "sim_qc-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args sce={input.sce} html={output}" {input.script} {log}'''

rule sim_data:
	priority: 98
	input:  script = config["scripts"] + "sim_data.R",
			sim_pars = config["sim_pars"] + "{sid}.json",
			sce = config["raw_data"] + "ref_{did}.rds"
	output: sim = config["sim_data"] + "{did},{sid},{i}.rds"
	log:	config["logs"] + "sim_data-{did},{sid},{i}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args sce={input.sce} sim={output.sim}\
		sim_pars={input.sim_pars} wcs={wildcards}"\
		{input.script} {log}'''

rule run_meth:
	priority: 97
	threads: 1
	input:	script = config["scripts"] + "run_meth.R",
			sim = config["sim_data"] + "{did},{sid},{i}.rds",
			meth_pars = config["meth_pars"] + "{mid}.json",
			run_pars = config["run_pars"] + "{did},{sid}.json",
			fun = lambda wc: config["scripts"] + "apply_" + mids.loc[wc.mid, "type"] + ".R"
	output: res = config["results"] + "{did},{sid},{i},{mid},{j},g{g},c{c},k{k},s{s}.rds"
	log:	config["logs"] + "run_meth-{did},{sid},{i},{mid},{j},g{g},c{c},k{k},s{s}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args sim={input.sim} fun={input.fun} wcs={wildcards}\
		meth_pars={input.meth_pars} run_pars={input.run_pars} res={output.res}"\
		{input.script} {log}'''

# pseudobulk-level mean-dispersion scatters
rule plot_pb_mean_disp:
	input:	config["utils"],
			script = config["scripts"] + "plot_pb_mean_disp.R",
			sce = config["raw_data"] + "ref_{did}.rds",
	output: ggp = config["plots"] + "{did}-pb_mean_disp.rds",
			fig = config["plots"] + "{did}-pb_mean_disp.pdf"
	log:	config["logs"] + "plot_pb_mean_disp-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args sce={input.sce} ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# null simulation p-value distributions
rule plot_null:
	input:	config["utils"],
			script = config["scripts"] + "plot_null.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",nill,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did}-null.rds", 
			fig = config["plots"] + "{did}-null.pdf"
	log:	config["logs"] + "plot_null-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR stratified by DD category using locally/globally adjusted p-values
rule plot_perf_by_cat:
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_cat.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z][0-9]+,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did}-perf_by_cat_{padj}.rds",
			fig = config["plots"] + "{did}-perf_by_cat_{padj}.pdf"
	log:	config["logs"] + "plot_perf_by_cat-{did},p_adj.{padj}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} wcs={wildcards}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR points stratified by #(cells)/#(samples)
rule plot_perf_by_nx:
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_nx.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10_n" + wc.x + ",.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did}-perf_by_n{x}.rds",
			fig = config["plots"] + "{did}-perf_by_n{x}.pdf"
	log:	config["logs"] + "plot_perf_by_n{x}-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} wcs={wildcards}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR points stratified by sample/group sizes
rule plot_perf_by_xs:
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_xs.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10_" + wc.y + "s,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did}-perf_by_{y}s.rds",
			fig = config["plots"] + "{did}-perf_by_{y}s.pdf"
	log:	config["logs"] + "plot_perf_by_{y}s-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR points stratified by expression
rule plot_perf_by_es:
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_es.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did}-perf_by_es_{padj}.rds",
			fig = config["plots"] + "{did}-perf_by_es_{padj}.pdf"
	log:	config["logs"] + "plot_perf_by_es-{did},p_adj.{padj}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} wcs={wildcards}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

rule plot_upset:
	input:	config["utils"],
			script = config["scripts"] + "plot_upset.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z][0-9]+,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: fig = config["plots"] + "{did}-upset.pdf"
	log:	config["logs"] + "plot_upset-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} fig={output.fig}"\
		{input.script} {log}'''

# scatter plots of estimated vs. simulated logFCs
rule plot_lfc:
	input:	config["utils"],
			script = config["scripts"] + "plot_lfc.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z]10,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output:	ggp = config["plots"] + "{did}-sim_vs_est_lfc.rds",
			fig = config["plots"] + "{did}-sim_vs_est_lfc.pdf"
	log:	config["logs"] + "plot_lfc-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# runtimes vs. #(cells)/#(genes)
rule plot_runtimes:
	input: 	config["utils"],
			script = config["scripts"] + "plot_runtimes.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10_n[g|c],*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output:	fig = config["plots"] + "{did}-runtimes.pdf"
	log:	config["logs"] + "plot_runtimes-{did}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} fig={output.fig}"\
		{input.script} {log}'''

#rule lps_run_meth:
#	input:	script = config["scripts"], "lps_run_meth.R"),
#			sce = config["raw_data"], "sce0_magl.rds"),
#			meth_pars = config["meth_pars"], "{mid}.json"),
#			fun = lambda wc: config["scripts"], "apply_" + mids.loc[wc.mid, "type"] + ".R")
#	output:	res = config["results"], "lps", "{mid}.rds")
#	script:	"{input.script}"

# write session info to .txt file
rule session_info:
	input:	config["scripts"] + "session_info.R"
	output:	"session_info.txt"
	log:	config["logs"] + "session_info.Rout" 
	shell:	'''{R} CMD BATCH --no-restore --no-save\
	"--args txt={output}" {input} {log}'''