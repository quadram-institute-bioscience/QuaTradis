import os.path
import sys
import subprocess
import shutil

from quatradis.util.parser import create_parser
import quatradis.util.file_handle_helpers as fhh


def add_subparser(subparsers):
    pipeline_parser_desc = "TraDIS pipelines that stitch together other tools in this package."
    pipeline_parser = subparsers.add_parser("pipeline", help=pipeline_parser_desc)
    pipeline_subparsers = pipeline_parser.add_subparsers(title=pipeline_parser_desc)

    create_parser("create_plots", pipeline_subparsers, create_plots_pipeline, create_plots_options,
                  "Creates transponson insertion site plot files for multiple fastqs in parallel where possible using snakemake.",
                  description='''This pipeline uses snakemake, therefore it is possible to customise how this operates to distribute the workload
    across a cluster using a snakemake config file.''',
                  usage="tradis pipeline create_plots [options] <fastq list file> <reference>")

    create_parser("compare", pipeline_subparsers, compare_pipeline, compare_options,
                  "Calculates gene essentiality for a set of transposon insertion site plots files ",
                  description='''This pipeline uses snakemake, therefore it is possible to customise how this operates to distribute the workload
    across a cluster using a snakemake config file.''',
                  usage="tradis pipeline compare [options]")


def create_plots_options(parser):
    parser.add_argument('fastqs', type=str,
                        help='Either a line separated text file containing a list of fastq formatted reads for processing or a single fastq file (can be gzipped).')
    parser.add_argument('reference', type=str,
                        help='The fasta formatted reference for processing.')
    parser.add_argument('-o', '--output_dir', default="results",
                        help='The output directory to use for all output files (default: results)')
    parser.add_argument('-n', '--threads', type=int, default=1,
                        help='number of threads to use when mapping and sorting (default: 1)')
    parser.add_argument('-a', '--aligner', default="bwa",
                        help='mapping tool to use (bwa, smalt, minimap2) (default: bwa)')
    parser.add_argument('-m', '--mapping_score', type=int, default=30,
                        help='mapping quality must be greater than X (Default: 30)')
    parser.add_argument('-t', '--tag', type=str, default="",
                        help='the tag to remove from fastq input (default: "")')
    parser.add_argument('-mm', '--mismatch', type=int, default=0,
                        help='number of mismatches allowed when matching tag (default: 0)')
    parser.add_argument('-sp', '--snakemake_profile', type=str,
                        help='If provided, pass this directory onto snakemake.  Assumes there is a file called "config.yaml" in that directory.')


def create_plots_pipeline(args):
    """
    Use snakemake to process multiple fastqs in parallel
    """

    fhh.ensure_output_dir_exists(args.output_dir)

    fastqs = [args.fastqs]
    if not fhh.is_fastq(args.fastqs):
        with open(args.fastqs, 'r') as fql:
            fastqs = [x.strip() for x in fql.readlines() if x]

    snakemake_config = os.path.join(args.output_dir, "create_plots_config.yaml")
    fastq_dir, fq_fn = os.path.split(args.fastqs)
    with open(snakemake_config, 'w') as ofql:
        ofql.write(create_yaml_option("output_dir", args.output_dir))
        ofql.write(create_yaml_option("reference", args.reference))
        ofql.write("fastq_dir: \"" + fastq_dir + "\"\n")
        ofql.write("fastqs:\n")
        for x in fastqs:
            ofql.write("- " + x + "\n")
        ofql.write(create_yaml_option("tag", args.tag))
        ofql.write(create_yaml_option("aligner", args.aligner))
        ofql.write(create_yaml_option("threads", args.threads, num=True))
        ofql.write(create_yaml_option("mismatch", args.mismatch, num=True))
        ofql.write(create_yaml_option("mapping_score", args.mapping_score, num=True))

    pipeline = find_pipeline_file("create_plots.smk")

    start_snakemake(pipeline, snakemake_config, threads=args.threads,
                                snakemake_profile=args.snakemake_profile, verbose=args.verbose)


def compare_options(parser):
    parser.add_argument('--condition_files', type=str, nargs='+',
                        help='A set of condition plot files to process (can be gzipped)')
    parser.add_argument('--control_files', type=str, nargs='+',
                        help='A set of control plot files to process (can be gzipped)')
    parser.add_argument('-o', '--output_dir', default="results",
                        help='The output directory to use for all output files (default: results)')
    parser.add_argument('-n', '--threads', type=int, default=1,
                        help='number of threads to use when processing (default: 1)')
    parser.add_argument('--annotations', '-a',
                        help='If provided genes in this EMBL annotations file will expanded based on data in the plotfile.',
                        type=str, default=None)
    parser.add_argument('--minimum_threshold', '-m',
                        help='Only include insert sites with this number or greater insertions', type=int, default=5)
    parser.add_argument('--prime_feature_size', '-z',
                        help='Feature size when adding 5/3 prime block when --use_annotation', type=int, default=198)
    parser.add_argument('--window_interval', '-l', help='Window interval', type=int, default=25)
    parser.add_argument('--window_size', '-w', help='Window size', type=int, default=100)
    parser.add_argument('-sp', '--snakemake_profile', type=str,
                        help='If provided, pass this directory onto snakemake.  Assumes there is a file called "config.yaml" in that directory.')


def compare_pipeline(args):
    """
    Use snakemake to process multiple fastqs in parallel
    """

    if len(args.control_files) != len(args.condition_files):
        raise ValueError("Must have equal number of control and condition files")

    if len(args.control_files) <= 1:
        raise ValueError("Must have 2 or more replicates of control and condition files")

    fhh.ensure_output_dir_exists(args.output_dir)

    snakemake_config = os.path.join(args.output_dir, "analysis_config.yaml")
    with open(snakemake_config, 'w') as ofql:
        ofql.write(create_yaml_option("output_dir", args.output_dir))
        ofql.write("condition_files:\n")
        for x in args.condition_files:
            ofql.write("- " + x + "\n")
        ofql.write("control_files:\n")
        for x in args.control_files:
            ofql.write("- " + x + "\n")
        ofql.write(create_yaml_option("threads", args.threads, num=True))
        ofql.write(create_yaml_option("annotations", args.annotations))
        ofql.write(create_yaml_option("minimum_threshold", args.minimum_threshold))
        ofql.write(create_yaml_option("prime_feature_size", args.prime_feature_size))
        ofql.write(create_yaml_option("window_interval", args.window_interval))
        ofql.write(create_yaml_option("window_size", args.window_size))

    pipeline = find_pipeline_file("compare.smk")

    start_snakemake(pipeline, snakemake_config, threads=args.threads, snakemake_profile=args.snakemake_profile,
                    verbose=args.verbose)


def start_snakemake(snakefile, snakemake_config, threads=1, snakemake_profile=None, verbose=False):
    if verbose:
        print("Using snakemake config files at: " + ", ".join([snakemake_config]))
        print("Starting snakemake pipeline")

    cmd_list = ["snakemake",
                "--snakefile=" + snakefile,
                "--configfile=" + snakemake_config,
                "--cores=" + str(threads),
                "--printshellcmds"]

    if snakemake_profile:
        cmd_list.append("--profile=" + snakemake_profile)
    if verbose:
        cmd_list.append("--verbose")

    cmd = " ".join(cmd_list)
    if verbose:
        print("Snakemake command:", cmd)

    exit_code = subprocess.call(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)

    if verbose:
        print("Snakemake returned exitcode:", exit_code)

    if exit_code != 0:
        raise RuntimeError("Error running pipeline.  Return code: " + str(exit_code))


def find_pipeline_file(pipeline_file):
    """
    Depending on how quatradis gets installed we may need to look in different locations for the snakemake pipeline file
    """
    local_path = os.path.join(os.path.dirname(__file__), pipeline_file)

    if os.path.exists(local_path):
        return local_path

    docker_path = os.path.join("/quatradis", "../../quatradis/pipelines", pipeline_file)

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
