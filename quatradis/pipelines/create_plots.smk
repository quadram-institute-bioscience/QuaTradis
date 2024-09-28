import os

seqs=[]

rule finish:
    input:
        stats=os.path.join(config["output_dir"], "combined.plot.stats"),
        refer_idx=os.path.join(config["output_dir"], os.path.basename(config["reference"]) + ".fai")
    run:
        print("All done!")


rule mapper_index_reference:
    input: config["reference"]
    output: touch(os.path.join(config["output_dir"],"index.done"))
    params:
        refname="myref"
    message: "Indexing reference"
    shell:
        "tradis utils index {input} {params.refname}"

rule samtools_index_reference:
    input: config["reference"]
    output: os.path.join(config["output_dir"], os.path.basename(config["reference"]) + ".fai")
    message: "Extracting sequence names from reference"
    shell:
        "samtools faidx {input} && cp {input}.fai {output}"

rule create_plot:
    input:
        idx=rules.mapper_index_reference.output,
        fq=os.path.join(config["fastq_dir"], "{fq}"),
        ref=config["reference"]
    output:
        stats=os.path.join(config["output_dir"], "{fq}", "tradis_out.plot.stats")
    params:
        aligner=("--aligner=" + config["aligner"]) if config["aligner"] else "",
        tag=("--tag=" + config["tag"]) if config["tag"] else "",
        mismatch=("--mismatch=" + str(config["mismatch"])) if config["mismatch"] else "",
        mapping_score="--mapping_score=" + str(config["mapping_score"]),
        threads="--threads=" + str(config["threads"]) if config["threads"] else "",
        output_dir=os.path.join(config["output_dir"], "{fq}")
    threads: int(config["threads"])
    message: "Creating transposon insertion site plot file for {input.fq}"
    shell:
        """tradis plot create {params.aligner} {params.threads} {params.tag} {params.mismatch} {params.mapping_score} \
        --output_dir={params.output_dir} --output_prefix=tradis_out --no_ref_index {input.fq} {input.ref}"""


rule combine_stats:
    input:
        stats=expand(os.path.join(config["output_dir"], "{fq}", "tradis_out.plot.stats"), fq=config["fastqs"])
    output: os.path.join(config["output_dir"], "combined.plot.stats")
    message: "Combining plot stats"
    shell:
        "head -n 1 {input[0]} > {output} && for i in {input}; do tail -n+2 $i >> {output} && echo "" >> {output}; done"

