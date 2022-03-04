import os.path
import sys
import time
import shutil

from quatradis.tags import remove_tags
from quatradis.mapper import calc_read_length, index_reference, map_reads, sam2bam
from quatradis.isp_create import plot


def run_tradis(fastq, reference, output_prefix, tag="", mapper="bwa", index=True, threads=1, max_mismatches=0, cutoff=30, verbose=False):
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
    :param max_mismatches: number of mismatches allowed when matching tag
    :param cutoff: Quality score cutoff value.  Alignments with score less than this are not considered for analysis.
    :param verbose: Extra logging information
    :return:
    """

    if not os.path.exists(fastq):
        raise ValueError("Fastq input file not found at " + fastq)

    if not os.path.exists(reference):
        raise ValueError("Reference fasta file not found at " + reference)

    start = time.time()

    outdir, prefix = os.path.split(output_prefix)

    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    detagged = output_prefix + ".rmtag.fastq.gz"

    if verbose:
        print("::::::::::::::::::\n" + fastq + "\n::::::::::::::::::\n\n", file=sys.stderr)

    nb_reads = 0
    nb_output_reads = 0
    nb_tagged_reads = 0
    if tag:
        if verbose:
            print("..........Removing tags that match user input: " + tag + "\n", file=sys.stderr)
        nb_reads, nb_output_reads, nb_tagged_reads = remove_tags(fastq, detagged, tag=tag, max_mismatches=max_mismatches, filter=True, trim=True)
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
    plot(bam, fastq, plot_file, cutoff_score=cutoff, nb_reads=nb_reads, nb_tagged_reads=nb_tagged_reads)

    end = time.time()
    if verbose:
        print("Tradis has completed in", '{:.3f}'.format(end - start) + "s")
        print("Mapped reads are here:", bam)
        print("Plot files are here:", plot_file + ".*")


def find_nextflow_file():
    """
    Depending on how quatradis gets installed we may need to look in different locations for the nextflow pipeline file
    """

    local_path = os.path.join(os.path.dirname(__file__), "..", "pipelines", "multi_tradis.nf")

    if os.path.exists(local_path):
        return local_path

    docker_path = os.path.join("/quatradis", "pipelines", "multi_tradis.nf")

    if os.path.exists(docker_path):
        return docker_path

    exe_path = shutil.which("multi_tradis.nf")

    if os.path.exists(exe_path):
        return exe_path

    raise RuntimeError("Could not find nextflow pipeline file.")


def run_multi_tradis(fastqs, reference, output_dir="results", nextflow_config="", tag="", aligner="bwa", threads=1, max_mismatches=0, cutoff=30, verbose=False):
    """
    Use nextflow to process multiple fastqs in parallel
    """

    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    with open(fastqs, 'r') as fql:
        fastqs = [x.strip() for x in fql.readlines() if x]
        cleaned_fastq_list = os.path.join(output_dir, "quadtradis_nf.fastq.txt")
        with open(cleaned_fastq_list, 'w') as ofql:
            ofql.write("\n".join(fastqs))

    nfcfg_cmd = ""
    if nextflow_config:
        nfcfg_cmd = "-c " + nextflow_config

    nf_pipeline = find_nextflow_file()

    command_args = ["nextflow", nfcfg_cmd, nf_pipeline, "--reference=" + reference,
               "--fastqs=" + cleaned_fastq_list, "--refname=myref", "--outdir=" + output_dir,
               ("--tag=" + tag) if tag else "", ("--aligner=" + aligner) if aligner else "",
               ("--threads=" + str(threads)) if threads else "",
               ("--mismatch=" + str(max_mismatches)) if max_mismatches else "",
               ("--mapping_score=" + str(cutoff)) if cutoff else ""]

    command = " ".join(command_args)

    if verbose:
        print("Attempting to execute nextflow pipeline with command:", command)

    os.system(command)
