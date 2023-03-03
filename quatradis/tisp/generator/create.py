import os.path
import sys
import time

import quatradis.tisp.generator.from_alignments as tisp_create

from quatradis.util.tags import remove_tags
from quatradis.util.mapper import calc_read_length, index_reference, map_reads, sam2bam
from quatradis.util import file_handle_helpers


def run_tradis(fastq, reference, output_prefix, alignments="", tag="", mapper="bwa", index=True, threads=1, mismatch=0, mapping_score=30, verbose=False):
    """
    The main program logic for running tradis.  Tradis takes in a fastq file and a reference as input and generates
    transposon insertion site plot files for each reference sequence, along with stats related to this.
    :param input: Either a fastq formatted file containing reads (can be read gzipped or uncompressed) or an alignment file (e.g. SAM/BAM).
    :param reference: The fasta formatted reference sequences
    :param output_prefix: Output prefix for output files.  If empty, then we output to current directory and use the fastq filename as prefix
    :param alignments: The alignment file (SAM or BAM) to use as input.  If provided, then we won't try to recreate the alignments directly from fastq.
    :param tag: Tag to identify and trim in reads.  If left empty then we run tradis in tagless mode and map reads as is.
    :param mapper: The mapping tools to map reads to reference (bwa, smalt, minimap2, minimap2_long)
    :param index: Whether to index the reference.  Set to false if reference is already indexed using the specified mapper.
    :param threads: Number of threads used for mapping and sorting alignments
    :param mismatch: number of mismatches allowed when matching tag
    :param mapping_score: Quality score cutoff value.  Alignments with score less than this are not considered for analysis.
    :param verbose: Extra logging information
    :return:
    """

    if not os.path.exists(fastq):
        raise ValueError("Input file not found at " + fastq)

    if not os.path.exists(reference):
        raise ValueError("Reference fasta file not found at " + reference)

    start = time.time()

    if not output_prefix:
        output_prefix = os.path.join(".", fastq)

    file_handle_helpers.ensure_output_dir_exists(output_prefix, includes_filename=True)

    alignment_file = ""
    nb_reads = 0
    nb_tagged_reads = 0
    if alignments:
        if not os.path.exists(alignments):
            raise ValueError("Alignments file not found at " + alignments)
        else:
            alignment_file = alignments
    else:
        alignment_file, nb_reads, nb_tagged_reads = create_alignments(fastq, reference, output_prefix, tag, mapper, index, threads, mismatch, verbose)


    plot_file = output_prefix + ".plot"
    if verbose:
        print("..........Generate insertion site plot files and statistics\n", file=sys.stderr)
    tisp_create.plot(alignment_file, fastq, plot_file, cutoff_score=mapping_score, nb_reads=nb_reads, nb_tagged_reads=nb_tagged_reads)

    end = time.time()
    if verbose:
        print("Tradis has completed in", '{:.3f}'.format(end - start) + "s")
        print("Mapped reads are here:", alignment_file)
        print("Plot files are here:", plot_file + ".*")



def create_alignments(input, reference, output_prefix, tag="", mapper="bwa", index=True, threads=1, mismatch=0, verbose=False):
    detagged = output_prefix + ".rmtag.fastq.gz"

    if verbose:
        print("::::::::::::::::::\n" + input + "\n::::::::::::::::::\n\n", file=sys.stderr)

    nb_reads = 0
    nb_tagged_reads = 0
    if tag:
        if verbose:
            print("..........Removing tags that match user input: " + tag + "\n", file=sys.stderr)
        nb_reads, nb_output_reads, nb_tagged_reads = remove_tags(input, detagged, tag=tag, max_mismatches=mismatch, filter=True, trim=True)
    else:
        if verbose:
            print("..........Tagless mode selected, skipping read preparation step\n", file=sys.stderr)
        detagged = input

    index_file = output_prefix + "ref.index"
    mapped_reads = output_prefix + ".mapped.sam"
    if verbose:
        print("..........Map reads to reference using " + mapper + " (using " + str(threads) + " threads)\n", file=sys.stderr)
    read_length = calc_read_length(detagged)
    if index:
        index_reference(reference, index_file, read_length, mapper)
    map_reads(detagged, reference, index_file, mapped_reads, read_length, mapper, threads)

    bam = output_prefix + ".mapped.bam"
    if verbose:
        print("..........Convert BAM, sort, index and check (using " + str(threads) + " threads)\n", file=sys.stderr)
    sam2bam(mapped_reads, bam, threads=threads)
    os.remove(mapped_reads)

    return bam, nb_reads, nb_tagged_reads