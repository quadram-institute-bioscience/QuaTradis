import os
import shutil
import sys

input_files=[]
norm_files=[]
plotnames=[]
controlnames=[]
conditionnames=[]
norm_lut={}
for p in config["condition_files"]:
    input_files.append(p)
    plotname = os.path.basename(p).split('.')[0]
    plotnames.append(plotname)
    conditionnames.append(plotname)


for p in config["control_files"]:
    input_files.append(p)
    plotname = os.path.basename(p).split('.')[0]
    plotnames.append(plotname)
    controlnames.append(plotname)

for p in input_files:
    plotname = os.path.basename(p).split('.')[0]
    n=os.path.join(config["output_dir"], "analysis", plotname, "original.plot.gz")
    norm_files.append(n)
    norm_lut[plotname]=n

ZCAT_CMD='gzcat'
if not shutil.which('gzcat'):
	if shutil.which('zcat'):
		ZCAT_CMD='zcat'
	else:
		raise Error("Couldn't find gzcat or zcat on your system.  Please install and try again.")

ZCAT_JOIN=' <' if sys.platform == "darwin" else ''

def make_no_normalise_cmd():
	cmds = []
	for i in range(len(input_files)):
		out_dir = os.path.join(config["output_dir"], "analysis", plotname)
		cmds.append("mkdir -p " + out_dir)
		cmds.append("if [[ $(file " + input_files[i] + " | grep ASCII | wc -l | cut -f 1) > 0 ]]; then gzip -c " + input_files[i] + " > " + norm_files[i] + "; else cp " + input_files[i] + " " + norm_files[i] + "; fi")

	return "; ".join(cmds)

def make_normalise_cmd():
    return "tradis plot normalise -o " + os.path.join(config["output_dir"], "analysis") + " -n original.plot.gz " + \
        ("--minimum_proportion_insertions=" + config["minimum_proportion_insertions"] if config["minimum_proportion_insertions"] else "") + \
            " " + " ".join(input_files)
            
# Print command function for easier traceability
def print_command(stage, input_files, output_files, params):
    print(f"Executing {stage} with:")
    print(f"  Inputs: {input_files}")
    print(f"  Outputs: {output_files}")
    print(f"  Params: {params}")

rule finish:
    input:
        expand(os.path.join(config["output_dir"], "gene_report.tsv")),
        expand(os.path.join(config["output_dir"], "comparison", "{type}", "plot_absscatter.png"), type=["original", "forward", "reverse", "combined"]),
        #expand(os.path.join(config["output_dir"], "comparison", "{type}", "{type}.compare.csv"), type=["original", "forward", "reverse", "combined"]),
        expand(os.path.join(config["output_dir"], "comparison", "{type}", "essentiality.csv"), type=["original", "forward", "reverse", "combined"]),
        expand(os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv.essen.csv"), type=["original", "combined", "forward", "reverse"], plot=plotnames)
    run:
        print("All done!")

# Modification 2.0
rule prepare_embl:
    input:
        plot=config["condition_files"][0],
        embl=config["annotations"]
    output:
        os.path.join(config["output_dir"], "prepared.embl")
    message:
        "Preparing embl annotations file"
    params:
        minimum_threshold="--minimum_threshold=" + str(config["minimum_threshold"]) if config.get("minimum_threshold") else "",
        window_size="--window_size=" + str(config["window_size"]) if config.get("window_size") else "",
        window_interval="--window_interval=" + str(config["window_interval"]) if config.get("window_interval") else "",
        prime_feature_size="--prime_feature_size=" + str(config["prime_feature_size"]) if config.get("prime_feature_size") else "",
        dynamic_window="--dynamic_window" if config.get("dynamic_window", False) else ""
    shell:
        "tradis compare prepare_embl --output={output} {params.minimum_threshold} {params.window_size} {params.window_interval} {params.prime_feature_size} {params.dynamic_window} --emblfile {input.embl} {input.plot}"

#    run:
#        print_command("prepare_embl", input, output, params)

rule normalise:
    input:
        input_files
    output:
        norm_files
    message: "Normalising plot files"
    log: os.path.join(config["output_dir"], "analysis", "normalise.log")
    params:
        cmd=make_normalise_cmd() if config["normalise"] else make_no_normalise_cmd()
    shell: "{params.cmd}"
#    run:
#        print_command("normalise", input, output, params)

rule split_plots:
    input:
        p=lambda wildcards: norm_lut[wildcards.plot]
    output:
        c=os.path.join(config["output_dir"], "analysis", "{plot}", "combined.plot.gz"),
        f=os.path.join(config["output_dir"], "analysis", "{plot}", "forward.plot.gz"),
        r=os.path.join(config["output_dir"], "analysis", "{plot}", "reverse.plot.gz")
    message: "Splitting plot file {input.p}"
    params:
        output_dir=os.path.join(config["output_dir"], "analysis", "{plot}"),
        minimum_threshold="--minimum_threshold=" + config["minimum_threshold"] if config["minimum_threshold"] else ""
    shell:
        """
        cmd="tradis compare split -o {params.output_dir} --gzipped {params.minimum_threshold} {input.p}"
        echo "Executing: $cmd"
        eval $cmd
        """
#    run:
#        print_command("split_plots", input, output, params)


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
    shell:
        """
        cmd="tradis plot count -o {params.output_dir} -s {params.suffix} {input.embl} {input.p}"
        echo "Executing: $cmd"
        eval $cmd
        """
#    run:
#        print_command("count_plots", input, output, params)


rule essentiality:
    input:
        c=os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv"),
    output:
        ess=os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv.essen.csv")
    log: os.path.join(config["output_dir"], "analysis", "{plot}", "{type}.count.tsv.essen.log")
    message: "Determining gene essentiality for {input}"
    shell:
        """
        cmd="tradis compare essentiality --verbose {input} > {log} 2>&1"
        echo "Executing: $cmd"
        eval $cmd
        """
#    run:
#        print_command("essentiality", input, output, {})


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
    shell:
        """
        cmd="tradis compare figures {params.window_size} {params.prefix} {params.controls} {params.conditions} > /dev/null 2>&1"
        echo "Executing: $cmd"
        eval $cmd
        """


rule insertion_site_comparison:
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
        span_gaps="--span_gaps=" + config["span_gaps"],
        minimum_block="--minimum_block=" + config["minimum_block"],
        minimum_logfc="--minimum_logfc=" + config["minimum_logfc"],
        minimum_logcpm="--minimum_logcpm=" + config["minimum_logcpm"],
        minimum_proportion_insertions="--minimum_proportion_insertions=" + config["minimum_proportion_insertions"],
        p_value="--pvalue=" + config["p_value"],
        q_value="--qvalue=" + config["q_value"],
        window_size="--window_size=" + config["window_size"],
        zcat=ZCAT_CMD + ZCAT_JOIN
    message: "Calculating logfc for {output}"
    shell:
        """
        SEQLENGTH=$({params.zcat} {params.combined} | wc -l | awk '{{ $1=$1 }};1')
        cmd="tradis compare insertion_sites {input.embl} {params.tt} {params.output_dir} --genome_length=${{SEQLENGTH}} {params.span_gaps} {params.window_size} {params.p_value} {params.q_value} {params.minimum_block} {params.minimum_logfc} {params.minimum_logcpm} {params.minimum_proportion_insertions} --verbose --controls {input.controls} --conditions {input.conditions} > {log} 2>&1"
        echo "Executing: $cmd"
        eval $cmd
        """


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
    shell:
        """
        cmd="tradis compare gene_report --combined={input.combined} --forward={input.forward} --reverse={input.rev} {params.scores} {params.window_size} {params.output_dir} {params.annotations} {input.embl}"
        echo "Executing: $cmd"
        eval $cmd
        """


rule essentiality_analysis:
    input:
        controls = lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{control}",wildcards.type + ".count.tsv"),control=controlnames),
        conditions = lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{condition}",wildcards.type + ".count.tsv"),condition=conditionnames),
        ess_controls= lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{control}",wildcards.type + ".count.tsv.essen.csv"),control=controlnames),
        ess_conditions=lambda wildcards: expand(os.path.join(config["output_dir"],"analysis","{condition}",wildcards.type + ".count.tsv.essen.csv"),condition=conditionnames)
    output:
        os.path.join(config["output_dir"], "comparison", "{type}", "essentiality.csv")
    log: os.path.join(config["output_dir"],"comparison","{type}","essentiality_analysis.log"),
    params:
        output_dir = "--output_dir=" + os.path.join(config["output_dir"], "comparison", "{type}"),
        t = "--type={type}"
    message: "Running essentiality analysis for {output}"
    shell:
        """
        cmd="tradis compare essentiality_analysis {params.output_dir} {params.t} --verbose --controls {input.controls} --conditions {input.conditions} --ess_controls {input.ess_controls} --ess_conditions {input.ess_conditions} > {log} 2>&1"
        echo "Executing: $cmd"
        eval $cmd
        """