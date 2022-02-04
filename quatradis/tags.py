"""
Preparing reads
"""
import os
import re

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from pysam.libcalignedsegment import CIGAR_OPS

from quatradis import file_handle_helpers


def add_tags(input_alignments, output_alignments=""):
    """
    Takes a SAM/BAM/CRAM file and creates a new BAM with tr and tq tags added to the sequence and quality strings.
    """
    if not output_alignments:
        output_alignments = os.path.splitext(input_alignments)[0] + ".tr.bam"

    with pysam.AlignmentFile(input_alignments, mode=file_handle_helpers.input_alignment_mode(input_alignments)) as alignments_in:
        with pysam.AlignmentFile(output_alignments, mode=file_handle_helpers.output_alignment_mode(output_alignments), header=alignments_in.header) as alignments_out:
            for a_in in alignments_in.fetch():
                a_tagged = add_tags_to_alignment(a_in)
                alignments_out.write(a_tagged)


def add_tags_to_alignment(alignment):
    new_alignment = alignment.__copy__()

    tr_tag = ""
    tq_tag = ""
    for t in alignment.tags:
        if t[0] == "tr":
            tr_tag = t[1]
        elif t[0] == "tq":
            tq_tag = t[1]

    if (not alignment.is_unmapped) and alignment.is_reverse:
        new_alignment.query_sequence = str(Seq(tr_tag + alignment.get_forward_sequence()).reverse_complement())
        new_alignment.query_qualities = pysam.qualitystring_to_array(tq_tag + pysam.array_to_qualitystring(alignment.get_forward_qualities()))[::-1] # [::-1] reverses the string
    else:
        new_alignment.query_sequence = tr_tag + alignment.query_sequence
        new_alignment.query_qualities = pysam.qualitystring_to_array(tq_tag + pysam.array_to_qualitystring(alignment.query_qualities))

    if not alignment.is_unmapped:
        new_alignment.cigartuples = [(CIGAR_OPS.CMATCH, new_alignment.query_length)]

    return new_alignment


def remove_tags(seq_file_in, seq_file_out="", tag="", max_mismatches=0, filter=True, trim=True):
    """
    Reads in a fastq file with tradis tags already attached to the start of the sequence.  Input file can be gzipped.
    Optionally removes reads that do not contain the provided tag at the start of the sequence.  Also optionally removes
    the tag from each read.  Default is to do both.
    Note that this logic is not the same as standard adapter trimming logic which would run a simple alignment process across the entire read and trim from there.
    Outputs a file *.tag.fastq.gz unless an alternative outfile name is specified.  Uncompressed fastq output is supported.
    :param seq_file_in: The filename / path to sequence file to be prepared
    :param seq_file_out: The filename / path to prepared sequence file.  Default is "<inputFileNamePrefix>.rmtag.fastq.gz>"
    :param tag: Any tag to search for
    :param max_mismatches: The maximum number of mismatches allowed when searching for the tag
    :param filter: Whether to discard reads that do no contain a tag (default: true)
    :param trim: Whether to trim tags off of reads (default: true)
    :return: Nothing
    """

    if seq_file_out == "":
        if seq_file_in.endswith('.gz'):
            seq_file_out = os.path.splitext(seq_file_in)[0]
        seq_file_out = os.path.splitext(seq_file_out)[0]
        seq_file_out += ".rmtag.fastq.gz"

    input_opener, input_mode = file_handle_helpers.reader_opener(seq_file_in)
    output_opener, output_mode = file_handle_helpers.writer_opener(seq_file_out)

    nb_input_reads = 0
    nb_output_reads = 0
    nb_tagged_reads = 0

    # Open the fastq file and process each record
    with input_opener(seq_file_in, input_mode) as input_handle:
        with output_opener(seq_file_out, output_mode) as output_handle:
            recs = SeqIO.parse(input_handle, "fastq")
            for rec in recs:
                nb_input_reads += 1
                if find_tag(str(rec.seq), tag, max_mismatches):
                    # If tag is found and we want to trim then we simply trim the start of the read
                    output_rec = rec[len(tag):] if trim else rec
                    SeqIO.write(output_rec, handle=output_handle, format='fastq')
                    nb_output_reads += 1
                    nb_tagged_reads += 1

                # If we are not filtering then keep this read with no tag in the output
                if not filter:
                    SeqIO.write(rec, handle=output_handle, format='fastq')
                    nb_output_reads += 1

    # Close file handles
    input_handle.close()
    if not seq_file_out.endswith('.gz'):
        output_handle.close()

    return nb_input_reads, nb_output_reads, nb_tagged_reads


def find_tag(seq, tag, max_mismatches=0):
    """
    Tries to identify the given tag within the given sequence, allowing for the given number of mismatches.
    We try to keep to the same logic as Bio-Tradis here.  That means if no mismatches allowed then
    we use a reg ex match.  This does however allow for use of wildcards etc.  If mismatches are
    allowed then we use our own logic.  We might want to revisit this at some point have have some
    logic that is close to typical adapter trimming, e.g. but running a local alignment.
    :param seq: The sequence to be analysed for the tag
    :param tag: The tag to be found
    :param max_mismatches: Number of mismatches tolerated
    :return: True if tag found, false otherwise
    """

    # First run some sanity checks to make sure tag exists and is of sufficient length and isn't longer than the seq,
    # otherwise return false
    if len(tag) < 3:
        return False
    if len(tag) > len(str(seq)):
        return False

    if max_mismatches == 0:
        match = re.match(tag, seq)
        return True if match else False
    else:
        mismatches = 0
        for i, c in enumerate(str(tag)):
            if c != seq[i]:
                mismatches += 1
            if mismatches > max_mismatches:
                return False
        return True


def tags_in_alignment(mapped_reads):
    """
    Simply checks the first alignment in a SAM/BAM/CRAM file for a "tr" tradis tag.  If present return True.
    :param mapped_reads: The SAM/BAM/CRAM file to check for a tradis tag
    :return: True is tradis tag is present in alignments, false otherwise
    """
    mode = file_handle_helpers.input_alignment_mode(mapped_reads)

    alignments = pysam.AlignmentFile(mapped_reads, mode)
    #TODO check this doesn't take ages for big files
    for read in alignments:
        for t in read.tags:
            if t[0] == "tr":
                return True
        break

    return False
