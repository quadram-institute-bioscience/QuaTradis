import os.path
import sys
import time
import shutil
import snakemake

from quatradis.tags import remove_tags
from quatradis.mapper import calc_read_length, index_reference, map_reads, sam2bam
from quatradis.isp_create import plot
from quatradis import file_handle_helpers


def run_tradis(fastq, reference, output_prefix, tag="", mapper="bwa", index=True, threads=1, mismatch=0, mapping_score=30, verbose=False):
    """
    The main program logic for running tradis.  Tradis takes in a fastq file and a reference as input and generates
    transposon insertion site plot files for each reference sequence, along with stats related to this.
    :param fastq: A fastq formatted file containing reads. Can be read gzipped or uncompressed.
    :param reference: The fasta formatted reference sequences
    :param output_prefix: Output prefix for output files
    :param tag: Tag to identify and trim in reads.  If left empty then we run tradis in tagless mode and map reads as is.
    :param mapper: The mapping tools to map reads to reference (bwa, smalt, minimap2, minimap2_long)
    :param index: Whether or not to index the reference.  Set to false if reference is already indexed using the specified mapper.
    :param threads: Number of threads used for mapping and sorting alignments
    :param mismatch: number of mismatches allowed when matching tag
    :param mapping_score: Quality score cutoff value.  Alignments with score less than this are not considered for analysis.
    :param verbose: Extra logging information
    :return:
    """

    if not os.path.exists(fastq):
        raise ValueError("Fastq input file not found at " + fastq)

    if not os.path.exists(reference):
        raise ValueError("Reference fasta file not found at " + reference)

    start = time.time()

    file_handle_helpers.ensure_output_dir_exists(output_prefix, includes_filename=True)

    detagged = output_prefix + ".rmtag.fastq.gz"

    if verbose:
        print("::::::::::::::::::\n" + fastq + "\n::::::::::::::::::\n\n", file=sys.stderr)

    nb_reads = 0
    nb_output_reads = 0
    nb_tagged_reads = 0
    if tag:
        if verbose:
            print("..........Removing tags that match user input: " + tag + "\n", file=sys.stderr)
        nb_reads, nb_output_reads, nb_tagged_reads = remove_tags(fastq, detagged, tag=tag, max_mismatches=mismatch, filter=True, trim=True)
    else:
        if verbose:
            print("..........Tagless mode selected, skipping read preparation step\n", file=sys.stderr)
        detagged = fastq

    index = output_prefix + "ref.index"
    mapped_reads = output_prefix + ".mapped.sam"
    if verbose:
        print("..........Map reads to reference using " + mapper + " (using " + str(threads) + " threads)\n", file=sys.stderr)
    read_length = calc_read_length(detagged)
    if index:
        index_reference(reference, index, read_length, mapper)
    map_reads(detagged, reference, index, mapped_reads, read_length, mapper, threads)

    bam = output_prefix + ".mapped.bam"
    if verbose:
        print("..........Convert BAM, sort, index and check (using " + str(threads) + " threads)\n", file=sys.stderr)
    sam2bam(mapped_reads, bam, threads=threads)
    os.remove(mapped_reads)

    plot_file = output_prefix + ".plot"
    if verbose:
        print("..........Generate insertion site plot files and statistics\n", file=sys.stderr)
    plot(bam, fastq, plot_file, cutoff_score=mapping_score, nb_reads=nb_reads, nb_tagged_reads=nb_tagged_reads)

    end = time.time()
    if verbose:
        print("Tradis has completed in", '{:.3f}'.format(end - start) + "s")
        print("Mapped reads are here:", bam)
        print("Plot files are here:", plot_file + ".*")


def find_pipeline_file():
    """
    Depending on how quatradis gets installed we may need to look in different locations for the nextflow pipeline file
    """
    pipeline_file = "Snakefile"
    local_path = os.path.join(os.path.dirname(__file__), "..", "pipelines", pipeline_file)

    if os.path.exists(local_path):
        return local_path

    docker_path = os.path.join("/quatradis", "pipelines", pipeline_file)

    if os.path.exists(docker_path):
        return docker_path

    exe_path = shutil.which(pipeline_file)

    if os.path.exists(exe_path):
        return exe_path

    raise RuntimeError("Could not find nextflow pipeline file.")


def create_yaml_option(option, value, num=False):
    opt = option + ": "
    if num:
        opt += str(value)
    else:
        opt += "\"" + str(value) + "\""
    opt += "\n"
    return opt


def run_multi_tradis(fastq_list, reference, output_dir="results", tag="", aligner="bwa", threads=1, mismatch=0, mapping_score=30, verbose=False,
                     snakemake_profile=None):
    """
    Use snakemake to process multiple fastqs in parallel
    """

    file_handle_helpers.ensure_output_dir_exists(output_dir)

    with open(fastq_list, 'r') as fql:
        fastqs = [x.strip() for x in fql.readlines() if x]
        snakemake_config = os.path.join(output_dir, "quadtradis_snakemake.yaml")
        with open(snakemake_config, 'w') as ofql:
            ofql.write(create_yaml_option("output_dir", output_dir))
            ofql.write(create_yaml_option("reference", reference))
            ofql.write("fastqs:\n")
            for x in fastqs:
                ofql.write("- " + x + "\n")
            ofql.write(create_yaml_option("tag", tag))
            ofql.write(create_yaml_option("aligner", aligner))
            ofql.write(create_yaml_option("threads", threads, num=True))
            ofql.write(create_yaml_option("mismatch", mismatch, num=True))
            ofql.write(create_yaml_option("mapping_score", mapping_score, num=True))

    pipeline = find_pipeline_file()

    config_files = [snakemake_config]


    if verbose:
        print("Using snakemake config files at: " + ", ".join(config_files))
        print("Starting snakemake pipeline")

    cmd_list = ["snakemake",
         "--snakefile=" + pipeline,
         "--configfile=" + snakemake_config,
         "--cores=" + str(threads)]

    if snakemake_profile:
        cmd_list.append("--profile=" + snakemake_profile)
    if verbose:
        cmd_list.append("--verbose")

    cmd = " ".join(cmd_list)
    if verbose:
        print("Snakemake command:", cmd)

    os.system(cmd)
