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
		# exclude all treat methods for now
		res_dirs.append(expand(\
			config["results"] + "{did},{sid},{i},{mid},{j},g{g},c{c},k{k},s{s}.rds",
			did = did, sid = sid, mid = mids.id[mids.id.str.find("treat") == -1],\
			i = range(1, sim_pars["nr"][0] + 1),\
			j = range(1, run_pars["nr"][0] + 1),\
			g = run_pars["ng"], c = run_pars["nc"],\
			k = run_pars["nk"], s = run_pars["ns"]))
		# include treat methods on for Kang as ref & simulation IDs dx10
		if bool(re.match(r"d[a-z][0-9]+$", sid)) and did == "kang":
			res_dirs.append(expand(\
				config["results"] + "{did},{sid},{i},{mid},{j},g{g},c{c},k{k},s{s}.rds",
				did = did, sid = sid, mid = mids.id[mids.id.str.find("treat") != -1],\
				i = range(1, sim_pars["nr"][0] + 1),\
				j = range(1, run_pars["nr"][0] + 1),\
				g = run_pars["ng"], c = run_pars["nc"],\
				k = run_pars["nk"], s = run_pars["ns"]))

sim_dirs = sum(sim_dirs, [])
res_dirs = sum(res_dirs, [])

subs = ["all", "treat"]

rule all: 
	input:	sim_dirs, res_dirs,
			#expand(config["raw_data"], "sce_{did}.rds"), did = config["dids"]),
			#expand(config["raw_data"], "ref_{did}.rds"), did = config["dids"]),
		# 'countsimQC' reports
			expand(config["plots"] + "{did},{inc}-sim_qc.html",\
				did = config["dids"], inc = subs),
		# null simulation p-value distributions 
			expand(config["plots"] + "{did},{inc}-null.{ext}",\
				did = config["dids"], inc = subs,\
				ext = ["rds", "pdf"]),
		# expand(config["plots"], "{did}_qc.html"), did = config["dids"]),
			expand(config["plots"] + "{did},{inc}-pb_mean_disp.pdf",\
				did = config["dids"], inc = subs),
			expand(config["plots"] + "{did},{inc}-upset.pdf",\
				did = config["dids"], inc = subs),
		# TPR-FDR stratified by DD category using locally/globally adjusted p-values
			expand(config["plots"] + "{did},{inc}-perf_by_cat_{padj}.{ext}",\
				did = config["dids"], inc = subs,\
				padj = ["loc", "glb"], ext = ["rds", "pdf"]),
		# scatters of simulated vs. estimated logFCs
			expand(config["plots"] + "{did},{inc}-sim_vs_est_lfc.{ext}",\
				did = config["dids"], inc = subs, ext = ["rds", "pdf"]),
		# TPR-FDR stratified by no. of cells per cluster-sample
			expand(config["plots"] + "{did},{inc}-perf_by_n{x}.{ext}",\
				did = config["dids"], inc = subs,\
				x = "c", ext = ["rds", "pdf"]),
		# TPR-FDR stratified by no. of replicates per group
			expand(config["plots"] + "{did},{inc}-perf_by_n{x}.{ext}",\
				did = "kang", inc = subs,\
				x = "s", ext = ["rds", "pdf"]),
		# TPR-FDR stratified by magnitude of sample-size unbalancing
			expand(config["plots"] + "{did},{inc}-perf_by_{x}s.{ext}",\
				did = config["dids"], inc = subs,\
				x = ["s", "g"], ext = ["rds", "pdf"]),
		# TPR-FDR stratified by expression level 
			expand(config["plots"] + "{did},{inc}-perf_by_es_{padj}.{ext}",\
				did = config["dids"], inc = subs,\
				padj = ["loc", "glb"], ext = ["rds", "pdf"]),
		# method runtimes versus no. cells/genes
			expand(config["plots"] + "{did},{inc}-runtimes.pdf", did = ["kang"])
			#expand(config["results"], "lps", "{mid}.rds"), mid = mids.id)
		# run all methods on LPS dataset
			#expand("LPS/output/DS_results_{mid}.rds", mid = mids.id)

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
	priority: 50
	input:	config["utils"],
			script = config["scripts"] + "plot_pb_mean_disp.R",
			sce = config["raw_data"] + "ref_{did}.rds",
	output: ggp = config["plots"] + "{did},{inc}-pb_mean_disp.rds",
			fig = config["plots"] + "{did},{inc}-pb_mean_disp.pdf"
	log:	config["logs"] + "plot_pb_mean_disp-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args sce={input.sce} ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# null simulation p-value distributions
rule plot_null:
	priority: 50
	input:	config["utils"],
			script = config["scripts"] + "plot_null.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",nill,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did},{inc}-null.rds", 
			fig = config["plots"] + "{did},{inc}-null.pdf"
	log:	config["logs"] + "plot_null-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR stratified by DD category using locally/globally adjusted p-values
rule plot_perf_by_cat:
	priority: 50
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_cat.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z][0-9]+,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did},{inc}-perf_by_cat_{padj}.rds",
			fig = config["plots"] + "{did},{inc}-perf_by_cat_{padj}.pdf"
	log:	config["logs"] + "plot_perf_by_cat-{did},{inc},p_adj.{padj}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} wcs={wildcards}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR points stratified by #(cells)/#(samples)
rule plot_perf_by_nx:
	priority: 50
	wildcard_constraints: x = "s|c"
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_nx.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10_n" + wc.x + ",.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did},{inc}-perf_by_n{x}.rds",
			fig = config["plots"] + "{did},{inc}-perf_by_n{x}.pdf"
	log:	config["logs"] + "plot_perf_by_n{x}-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} wcs={wildcards}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR points stratified by sample/group sizes
rule plot_perf_by_xs:
	priority: 50
	wildcard_constraints: x = "s|g"
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_xs.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10_" + wc.x + "s[0-9],.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did},{inc}-perf_by_{x}s.rds",
			fig = config["plots"] + "{did},{inc}-perf_by_{x}s.pdf"
	log:	config["logs"] + "plot_perf_by_{x}s-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# TPR-FDR points stratified by expression
rule plot_perf_by_es:
	priority: 50
	input:	config["utils"],
			script = config["scripts"] + "plot_perf_by_es.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z][0-9]+,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: ggp = config["plots"] + "{did},{inc}-perf_by_es_{padj}.rds",
			fig = config["plots"] + "{did},{inc}-perf_by_es_{padj}.pdf"
	log:	config["logs"] + "plot_perf_by_es-{did},{inc},p_adj.{padj}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} wcs={wildcards}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

rule plot_upset:
	priority: 50
	input:	config["utils"],
			script = config["scripts"] + "plot_upset.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z][0-9]+,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output: fig = config["plots"] + "{did},{inc}-upset.pdf"
	log:	config["logs"] + "plot_upset-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} fig={output.fig}"\
		{input.script} {log}'''

# scatter plots of estimated vs. simulated logFCs
rule plot_lfc:
	priority: 50
	input:	config["utils"],
			script = config["scripts"] + "plot_lfc.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",d[a-z]10,.*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output:	ggp = config["plots"] + "{did},{inc}-sim_vs_est_lfc.rds",
			fig = config["plots"] + "{did},{inc}-sim_vs_est_lfc.pdf"
	log:	config["logs"] + "plot_lfc-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res}\
		ggp={output.ggp} fig={output.fig}"\
		{input.script} {log}'''

# runtimes vs. #(cells)/#(genes)
rule plot_runtimes:
	priority: 50
	input: 	config["utils"],
			script = config["scripts"] + "plot_runtimes.R",
			res = lambda wc: filter(re.compile(\
				config["results"] + wc.did +\
				",de10_n[g|c],*").search, res_dirs)
	params:	res = lambda wc, input: ";".join(input.res)
	output:	fig = config["plots"] + "{did},{inc}-runtimes.pdf"
	log:	config["logs"] + "plot_runtimes-{did},{inc}.Rout"
	shell:	'''{R} CMD BATCH --no-restore --no-save\
		"--args res={params.res} fig={output.fig}"\
		{input.script} {log}'''

# rule run_meth_lps:
# 	threads: 20
# 	priority: 10
# 	wildcard_constraints: mid = ".+(?!treat).+"
# 	input:	script = config["scripts"] + "run_meth_lps.R",
# 			sce = "LPS/output/SCE_annotation.rds",
# 			meth_pars = config["meth_pars"] + "{mid}.json",
# 			fun = lambda wc: config["scripts"] + "apply_" + mids.loc[wc.mid, "type"] + ".R"
# 	output:	res = "LPS/output/DS_results_{mid}.rds"
# 	log:	config["logs"] + "run_meth_lps-{mid}.Rout"
# 	shell:	'''{R} CMD BATCH --no-restore --no-save\
# 		"--args sce={input.sce} fun={input.fun} wcs={wildcards}\
# 		meth_pars={input.meth_pars} res={output.res} n_threads='20'"\
# 		{input.script} {log}'''

# write session info to .txt file
rule session_info:
	priority: 1
	input:	config["scripts"] + "session_info.R"
	output:	"session_info.txt"
	log:	config["logs"] + "session_info.Rout" 
	shell:	'''{R} CMD BATCH --no-restore --no-save\
	"--args txt={output}" {input} {log}'''