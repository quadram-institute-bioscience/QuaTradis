import os.path
import sys
import time

from quatradis.tags import remove_tags
from quatradis.mapper import index_and_align, sam2bam
from quatradis.isp_create import plot


def run_tradis(fastq, reference, output_prefix, tag="", mapper="bwa", threads=1, max_mismatches=0, cutoff=30, verbose=False):
    """
    The main program logic for running tradis.  Tradis takes in a fastq file and a reference as input and generates
    transposon insertion site plot files for each reference sequence, along with stats related to this.
    :param fastq: A fastq formatted file containing reads. Can be read gzipped or uncompressed.
    :param reference: The fasta formatted reference sequences
    :param output_prefix: Output prefix for output files
    :param tag: Tag to identify and trim in reads.  If left empty then we run tradis in tagless mode and map reads as is.
    :param mapper: The mapping tools to map reads to reference (bwa, smalt, minimap2, minimap2_long)
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
    nb_tagged_reads = 0
    if tag:
        if verbose:
            print("..........Removing tags that match user input: " + tag + "\n", file=sys.stderr)
        nb_reads, nb_tagged_reads = remove_tags(fastq, detagged, tag=tag, max_mismatches=max_mismatches, filter=True, trim=True)
    else:
        if verbose:
            print("..........Tagless mode selected, skipping read preparation step\n", file=sys.stderr)
        detagged = fastq

    index = output_prefix + "ref.index"
    mapped_reads = output_prefix + ".mapped.sam"
    if verbose:
        print("..........Map reads to reference using " + mapper + " (using " + str(threads) + " threads)\n", file=sys.stderr)
    index_and_align(detagged, reference, index, mapped_reads, mapper, threads=threads)

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
