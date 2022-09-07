import os.path
import sys
import time

import quatradis.tisp.generator.from_alignments as tisp_create

from quatradis.util.tags import remove_tags
from quatradis.util.mapper import calc_read_length, index_reference, map_reads, sam2bam
from quatradis.util import file_handle_helpers


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
    tisp_create.plot(bam, fastq, plot_file, cutoff_score=mapping_score, nb_reads=nb_reads, nb_tagged_reads=nb_tagged_reads)

    end = time.time()
    if verbose:
        print("Tradis has completed in", '{:.3f}'.format(end - start) + "s")
        print("Mapped reads are here:", bam)
        print("Plot files are here:", plot_file + ".*")



