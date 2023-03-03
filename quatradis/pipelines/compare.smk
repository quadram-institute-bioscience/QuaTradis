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
        expand(os.path.join(config["output_dir"], "gene_report.csv")),
        expand(os.path.join(config["output_dir"], "figures", "{type}", "plot_absscatter.png"), type=["forward", "reverse", "combined"])
    run:
        print("All done!")


rule prepareEMBL:
    input:
        plot=config["control_files"][0]
    output: os.path.join(config["output_dir"], "prepared.embl")
    message: "Preparing embl annotations file"
    params:
        minimum_threshold="--minimum_threshold=" + config["minimum_threshold"] if config["minimum_threshold"] else "",
        window_size="--window_size=" + config["window_size"] if config["window_size"] else "",
        window_interval="--window_interval=" + config["window_interval"] if config["window_interval"] else "",
        prime_feature_size="--prime_feature_size=" + config["prime_feature_size"] if config["prime_feature_size"] else "",
        emblfile="--emblfile=" + config["annotations"] if not config["annotations"]=="None" else ""
    shell: "tradis utils prepare_embl --output={output} {params.minimum_threshold} {params.window_size} {params.window_interval} {params.prime_feature_size} {params.emblfile} {input.plot}"

rule essentiality:
    input:
        p=lambda wildcards: plot_lut[wildcards.plot],
        embl=rules.prepareEMBL.output
    output:
        os.path.join(config["output_dir"], "essentiality", "{plot}", "combined.plot.gz"),
        os.path.join(config["output_dir"], "essentiality", "{plot}", "forward.plot.gz"),
        os.path.join(config["output_dir"], "essentiality", "{plot}", "reverse.plot.gz")
    message: "Determining gene essentiality for {input}"
    params:
        output_dir="--output_dir=" + os.path.join(config["output_dir"], "essentiality", "{plot}".split('.')[0]),
        minimum_threshold="--minimum_threshold=" + config["minimum_threshold"] if config["minimum_threshold"] else ""
    shell: "tradis utils essentiality {params.output_dir} {params.minimum_threshold} {input.p} {input.embl}"


rule plot:
    input:
        controls=lambda wildcards: expand(os.path.join(config["output_dir"],"essentiality","{control}",wildcards.type + ".plot.gz"),control=controlnames),
        conditions=lambda wildcards: expand(os.path.join(config["output_dir"],"essentiality","{condition}",wildcards.type + ".plot.gz"),condition=conditionnames),
    output: os.path.join(config["output_dir"], "figures", "{type}", "plot_absscatter.png")
    params:
        typestr=lambda wildcards: wildcards.type,
        prefix=lambda wildcards: "--prefix=" + os.path.join(config["output_dir"], "figures", wildcards.type, "plot"),
        window_size="--window_size=" + config["window_size"],
        controls=lambda wildcards: "--controls " + " ".join(expand(os.path.join(config["output_dir"], "essentiality", "{control}", wildcards.type + ".plot.gz"), control=controlnames)),
        conditions=lambda wildcards: "--conditions " + " ".join(expand(os.path.join(config["output_dir"], "essentiality", "{condition}", wildcards.type + ".plot.gz"), condition=conditionnames))
    message: "Creating figures for type: {params.typestr}"
    shell: "tradis compare figures {params.window_size} {params.prefix} {params.controls} {params.conditions}"


rule logfc:
    input:
        controls=lambda wildcards: expand(os.path.join(config["output_dir"], "essentiality", "{control}", wildcards.type + ".plot.gz"), control=controlnames),
        conditions=lambda wildcards: expand(os.path.join(config["output_dir"], "essentiality", "{condition}", wildcards.type + ".plot.gz"), condition=conditionnames),
        embl=rules.prepareEMBL.output
    output:
        os.path.join(config["output_dir"], "logfc", "{type}.csv")
    params:
        type=lambda wildcards: wildcards.type,
        combined=os.path.join(config["output_dir"], "essentiality", controlnames[0], "combined.plot.gz"),
        output_dir="--prefix=" + os.path.join(config["output_dir"], "logfc"),
        control_dirs="--control_dirs " + " ".join(expand(os.path.join(config["output_dir"], "essentiality", "{control}"), control=controlnames)),
        condition_dirs="--condition_dirs " + " ".join(expand(os.path.join(config["output_dir"],"essentiality","{condition}"),condition=conditionnames)),
        zcat=ZCAT_CMD
    message: "Calculating logfc for {output}"
    shell: "SEQLENGTH=$({params.zcat} {params.combined} | wc -l); tradis compare logfc_plot {input.embl} {params.type} {params.output_dir} -g $SEQLENGTH {params.control_dirs} {params.condition_dirs}"


rule gene_stats:
    input:
        combined=expand(os.path.join(config["output_dir"], "logfc", "{type}.csv"), type=["forward", "reverse", "combined"]),
        embl=rules.prepareEMBL.output
    output: os.path.join(config["output_dir"], "gene_report.csv")
    params:
        input_dir=os.path.join(config["output_dir"], "logfc"),
        output_dir="--output_dir=" + config["output_dir"],
        window_size="--window_size=" + config["window_size"],
        annotations="--annotations=" + config["annotations"] if not config["annotations"]=="None" else ""
    message: "Creating gene report"
    shell: "tradis compare gene_report {params.window_size} {params.output_dir} {params.input_dir} {params.annotations} {input.embl}"