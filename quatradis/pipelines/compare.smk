import os
import shutil

plotnames=[]
controlnames=[]
conditionnames=[]
plot_lut={}
for p in config["condition_files"]:
    plotname = os.path.basename(p).split('.')[0]
    plotnames.append(plotname)
    conditionnames.append(plotname)
    plot_lut[plotname]=p

for p in config["control_files"]:
    plotname = os.path.basename(p).split('.')[0]
    plotnames.append(plotname)
    controlnames.append(plotname)
    plot_lut[plotname]=p

ZCAT_CMD='gzcat'
if not shutil.which('gzcat'):
	if shutil.which('zcat'):
		ZCAT_CMD='zcat'
	else:
		raise Error("Couldn't find gzcat or zcat on your system.  Please install and try again.")


rule finish:
    input:
        expand(os.path.join(config["output_dir"], "gene_report.tsv")),
        expand(os.path.join(config["output_dir"], "comparison", "{type}", "plot_absscatter.png"), type=["original", "forward", "reverse", "combined"]),
        #expand(os.path.join(config["output_dir"],"comparison","{type}","{type}.compare.csv"), type=["original", "forward", "reverse", "combined"]),
        expand(os.path.join(config["output_dir"],"comparison","{type}","essentiality.csv"), type=["original", "forward", "reverse", "combined"]),
        expand(os.path.join(config["output_dir"],"analysis","{plot}","{type}.count.tsv.all.csv"), type=["original", "combined", "forward", "reverse"], plot=plotnames)
    run:
        print("All done!")


rule prepare_embl:
    input:
        plot=config["control_files"][0],
        embl=config["annotations"]
    output: os.path.join(config["output_dir"], "prepared.embl")
    message: "Preparing embl annotations file"
    params:
        minimum_threshold="--minimum_threshold=" + config["minimum_threshold"] if config["minimum_threshold"] else "",
        window_size="--window_size=" + config["window_size"] if config["window_size"] else "",
        window_interval="--window_interval=" + config["window_interval"] if config["window_interval"] else "",
        prime_feature_size="--prime_feature_size=" + config["prime_feature_size"] if config["prime_feature_size"] else ""
    shell: "tradis compare prepare_embl --output={output} {params.minimum_threshold} {params.window_size} {params.window_interval} {params.prime_feature_size} --emblfile {input.embl} {input.plot}"


rule split_plots:
    input:
        p=lambda wildcards: plot_lut[wildcards.plot]
    output:
        c=os.path.join(config["output_dir"], "analysis", "{plot}", "combined.plot.gz"),
        f=os.path.join(config["output_dir"], "analysis", "{plot}", "forward.plot.gz"),
        r=os.path.join(config["output_dir"], "analysis", "{plot}", "reverse.plot.gz")
    message: "Splitting plot file {input.p}"
    params:
        output_dir=os.path.join(config["output_dir"], "analysis", "{plot}"),
        minimum_threshold="--minimum_threshold=" + config["minimum_threshold"] if config["minimum_threshold"] else ""
    shell: "tradis compare split -o {params.output_dir} --gzipped {params.minimum_threshold} {input.p}"


rule count_plots:
    input:
        p=os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.plot.gz"),
        embl=rules.prepare_embl.output
    output:
        os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv")
    message: "Analysing combined plot and prepared embl file {input.p}, {input.embl}"
    params:
        output_dir=os.path.join(config["output_dir"], "analysis", "{plot}"),
        suffix="count.tsv"
    wildcard_constraints:
        type="(?!original).*"
    shell: "tradis plot count -o {params.output_dir} -s {params.suffix} {input.embl} {input.p}"


rule copy_original:
    input:
        lambda wildcards: plot_lut[wildcards.plot]
    output:
        os.path.join(config["output_dir"],"analysis","{plot}","original.plot.gz")
    message: "Copying original plot file to target folder: {input}"
    params:
        output_dir=os.path.join(config["output_dir"],"analysis","{plot}")
    shell: "if [[ $(file {input} | grep ASCII | wc -l | cut -f 1) > 0 ]]; then gzip -c {input} > {output}; else cp {input} {output}; fi"


rule count_original_plots:
	input:
		p=rules.copy_original.output,
		embl=config["annotations"]
	output:
		os.path.join(config["output_dir"], "analysis", "{plot}", "original.count.tsv")
	message: "Analysing original plot and embl file {input.p}, {input.embl}"
	params:
		output_dir=os.path.join(config["output_dir"], "analysis", "{plot}")
	shell: "tradis plot count -o {params.output_dir} -s count.tsv {input.embl} {input.p}"


rule essentiality:
    input:
        c=os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv"),
    output:
        all=os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv.all.csv"),
        ess=os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv.essen.csv")
    log: os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv.essen.log")
    message: "Determining gene essentiality for {input}"
    shell: "tradis compare essentiality --verbose {input} > {log} 2>&1"


rule plot:
    input:
        controls=lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{control}",wildcards.type + ".plot.gz"),control=controlnames),
        conditions=lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{condition}",wildcards.type + ".plot.gz"),condition=conditionnames),
    output: os.path.join(config["output_dir"], "comparison", "{type}", "plot_absscatter.png")
    params:
        typestr=lambda wildcards: wildcards.type,
        prefix=lambda wildcards: "--prefix=" + os.path.join(config["output_dir"], "comparison", wildcards.type, "plot"),
        window_size="--window_size=" + config["window_size"],
        controls=lambda wildcards: "--controls " + " ".join(expand(os.path.join(config["output_dir"], "analysis", "{control}", wildcards.type + ".plot.gz"), control=controlnames)),
        conditions=lambda wildcards: "--conditions " + " ".join(expand(os.path.join(config["output_dir"], "analysis", "{condition}", wildcards.type + ".plot.gz"), condition=conditionnames))
    message: "Creating figures for type: {params.typestr}"
    shell: "tradis compare figures {params.window_size} {params.prefix} {params.controls} {params.conditions} > /dev/null 2>&1"


rule logfc:
    input:
        controls=lambda wildcards: expand(os.path.join(config["output_dir"], "analysis", "{control}", wildcards.tt + ".count.tsv"), control=controlnames),
        conditions=lambda wildcards: expand(os.path.join(config["output_dir"], "analysis", "{condition}", wildcards.tt + ".count.tsv"), condition=conditionnames),
        embl=rules.prepare_embl.output
    output:
        os.path.join(config["output_dir"], "comparison", "{tt}", "{tt}.compare.csv"),
        os.path.join(config["output_dir"],"comparison","{tt}","{tt}.logfc.plot"),
    log: os.path.join(config["output_dir"], "comparison", "{tt}", "logfc.log"),
    params:
        tt=lambda wildcards: wildcards.tt,
        combined=os.path.join(config["output_dir"], "analysis", controlnames[0], "combined.plot.gz"),
        output_dir="--prefix=" + os.path.join(config["output_dir"], "comparison", "{tt}"),
        zcat=ZCAT_CMD
    message: "Calculating logfc for {output}"
    shell: "SEQLENGTH=$({params.zcat} {params.combined} | wc -l); tradis compare logfc_plot {input.embl} {params.tt} {params.output_dir} -g ${{SEQLENGTH}} --verbose --controls {input.controls} --conditions {input.conditions} > {log} 2>&1"


rule gene_stats:
    input:
        combined=os.path.join(config["output_dir"], "comparison", "combined", "combined.logfc.plot"),
        forward=os.path.join(config["output_dir"],"comparison","forward","forward.logfc.plot"),
        rev=os.path.join(config["output_dir"],"comparison","reverse","reverse.logfc.plot"),
        embl=rules.prepare_embl.output
    output: os.path.join(config["output_dir"], "gene_report.tsv")
    params:
        input_dir=os.path.join(config["output_dir"], "comparison"),
        output_dir="--output_dir=" + config["output_dir"],
        window_size="--window_size=" + config["window_size"],
        annotations="--annotations=" + config["annotations"] if not config["annotations"]=="None" else "",
        scores="--scores=" + os.path.join(config["output_dir"],"comparison","combined","combined.pqvals.plot")
    message: "Creating gene report"
    shell: "tradis compare gene_report --combined={input.combined} --forward={input.forward} --reverse={input.rev} {params.scores} {params.window_size} {params.output_dir} {params.annotations} {input.embl}"


rule essentiality_analysis:
    input:
        controls = lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{control}",wildcards.type + ".count.tsv.all.csv"),control=controlnames),
        conditions = lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{condition}",wildcards.type + ".count.tsv.all.csv"),condition=conditionnames),
        ess_controls= lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{control}",wildcards.type + ".count.tsv.essen.csv"),control=controlnames),
        ess_conditions=lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{condition}",wildcards.type + ".count.tsv.essen.csv"),condition=conditionnames)
    output:
        os.path.join(config["output_dir"], "comparison", "{type}", "essentiality.csv")
    log: os.path.join(config["output_dir"],"comparison","{type}","essentiality_analysis.log"),
    params:
        output_dir = "--output_dir=" + os.path.join(config["output_dir"], "comparison", "{type}"),
        t = "--type={type}"
    message: "Running essentiality analysis for {output}"
    shell: "tradis compare essentiality_analysis {params.output_dir} {params.t} --verbose --controls {input.controls} --conditions {input.conditions} --ess_controls {input.ess_controls} --ess_conditions {input.ess_conditions} > {log} 2>&1"
