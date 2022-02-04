import gzip
import os

from Bio import bgzf


def reader_opener(fastq_file):
    """
    Creates an input file opener based on the given file's extension: gzip file opener if gzip compression is required, or normal
    open opener if uncompressed.  Also returns the correct write mode that is required for the opener.
    :param fastq_file: fastq file that we want to read from, can be gzip compressed if .gz extension is present
    :return: input opener, and opening mode
    """
    if fastq_file.endswith('.gz'):
        input_opener = gzip.open
        input_mode = "rt"
    else:
        input_opener = open
        input_mode = "r"

    return input_opener, input_mode


def writer_opener(fastq_file):
    """
    Creates an output file opener based on the given file's extension: bgz file opener if gzip compression is required, or normal
    open opener if uncompressed.  Also returns the correct write mode that is required for the opener.
    :param fastq_file: fastq file that we want to write to, can be gzip compressed if .gz extension is present
    :return: output opener, and opening mode
    """
    if fastq_file.endswith('.gz'):
        output_opener = bgzf.BgzfWriter
        output_mode = "wb"
    else:
        output_opener = open
        output_mode = "w"
    return output_opener, output_mode


def input_alignment_mode(alignment_file):
    """
    Opening SAM, BAM, or CRAM files requires a different opening mode.  This function returns the correct mode based on
    the file extension of the given alignment file.
    :param alignment_file: SAM / BAM /CRAM file as determined by its extension
    :return: opener mode suitable for reading from file type
    """
    ext = os.path.splitext(alignment_file)[1]
    if ext == '.bam':
        mode = "rb"
    elif ext == '.sam':
        mode = "r"
    elif ext == '.cram':
        mode = "rc"
    else:
        raise ValueError("Invalid alignment format: {}".format(ext))

    return mode


def output_alignment_mode(alignment_file):
    """
    Opening SAM, BAM, or CRAM files requires a different opening mode.  This function returns the correct mode based on
    the file extension of the given alignment file.
    :param alignment_file: SAM / BAM /CRAM file as determined by its extension
    :return: opener mode suitable for writing to this file type
    """
    ext = os.path.splitext(alignment_file)[1]
    if ext == '.bam':
        mode = "wb"
    elif ext == '.sam':
        mode = "w"
    elif ext == '.cram':
        mode = "wc"
    else:
        raise ValueError("Invalid alignment format: {}".format(ext))

    return mode
