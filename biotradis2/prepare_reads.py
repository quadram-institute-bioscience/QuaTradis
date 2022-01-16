#!/usr/bin/env python3
# coding: utf-8

"""
Preparing reads for biotradis2
"""


import os
from Bio import SeqIO, bgzf
import gzip


def prepare_reads(seq_file_in, seq_file_out="", tag="", max_mismatches=0):
    """
    Reads in a fastq file with tradis tags already attached to the start of the sequence.  Input file can be gzipped.
    Removes reads that do not contain the provided tag at the start of the sequence.  Note that this logic is different
    from biotradis 1 which would keep the reads if no tag was found.  Also note that this logic is not the same as
    standard adapter trimming logic which would run a simple alignment process across the entire read and trim from there.
    Outputs a file *.tag.fastq.gz unless an alternative outfile name is specified.  Uncompressed fastq output is supported.
    :param seq_file_in: The filename / path to sequence file to be prepared
    :param seq_file_out: The filename / path to prepared sequence file.  Default is "<inputFileNamePrefix>.rmtag.fastq.gz>"
    :param tag: Any tag to search for
    :param max_mismatches: The maximum number of mismatches allowed when searching for the tag
    :return: Nothing
    """

    if seq_file_out == "":
        if seq_file_in.endswith('.gz'):
            seq_file_out = os.path.splitext(seq_file_in)[0]
        seq_file_out = os.path.splitext(seq_file_out)[0]
        seq_file_out += ".rmtag.fastq.gz"

    # Make sure we can handle gzipped or non-gzipped fastq for both input
    if seq_file_in.endswith('.gz'):
        input_opener = gzip.open
        input_mode = "rt"
    else:
        input_opener = open
        input_mode = "r"

    # Make sure we can handle gzipped or non-gzipped fastq for both input
    if seq_file_out.endswith('.gz'):
        output_opener = bgzf.BgzfWriter
        output_mode = "wb"
    else:
        output_opener = open
        output_mode = "w"

    # Open the fastq file and process each record
    with input_opener(seq_file_in, input_mode) as inputHandle:
        with output_opener(seq_file_out, output_mode) as outputHandle:
            recs = SeqIO.parse(inputHandle, "fastq")
            for rec in recs:
                if find_tag(rec.seq, tag, max_mismatches):
                    # If tag is found we simply trim off the start of the read
                    trimmed = rec[len(tag):]
                    SeqIO.write(trimmed, handle=outputHandle, format='fastq')

    # Close file handles
    inputHandle.close()
    if not seq_file_out.endswith('.gz'):
        outputHandle.close()


def find_tag(seq, tag, max_mismatches=0):
    """
    Roughly following adaptor removal technique described here: https://bcbio.wordpress.com/2009/08/09/trimming-adaptors-from-short-read-sequences/
    :param seq: The sequence to be processed
    :param tag: The tag to be removed
    :param max_mismatches: Number of mismatches tolerated
    :return: Sequence with tag removed
    """

    # First run some sanity checks to make sure tag exists and is of sufficient length and isn't longer than the seq, otherwise return false
    if len(tag) < 3:
        return False
    if len(tag) > len(str(seq)):
        return False
    mismatches = 0
    for i, c in enumerate(str(tag)):
        if c != seq[i]:
            mismatches += 1
        if mismatches > max_mismatches:
            return False
    return True
