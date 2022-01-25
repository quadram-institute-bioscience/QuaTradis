"""
Mapping prepared reads to index
"""
import gzip
import os

import pysam
from Bio import SeqIO


def index_and_align(reads, reference, index, out_file, mapper="bwa", threads=1):
    read_length = calc_read_length(reads)
    index_cmd = index_reference(reference, index, read_length, mapper)
    align_cmd = map_reads(reads, reference, index, out_file, read_length, mapper, threads)
    return index_cmd, align_cmd

def calc_read_length(reads_file):
    # Make sure we can handle gzipped or non-gzipped fastq for both input
    if reads_file.endswith('.gz'):
        input_opener = gzip.open
        input_mode = "rt"
    else:
        input_opener = open
        input_mode = "r"

    read_length = 0
    with input_opener(reads_file, input_mode) as inputHandle:
        records = SeqIO.parse(inputHandle, "fastq")
        for rec in records:
            read_length = len(rec.seq)
            break

    inputHandle.close()

    return read_length


def index_reference(reference, refname, read_length, mapper="bwa", dry_run=False):
    if mapper == "smalt":
        k = smalt_k_default(read_length)
        s = smalt_s_default(read_length)
        index_cmd = "smalt index -k " + str(k) + " -s " + str(s) + " " + refname + " " + reference + " > /dev/null 2>&1"
    elif mapper == "bwa":
        index_cmd = "bwa index " + reference + " > /dev/null 2>&1"
    elif mapper == "minimap2" or mapper == "minimap2_long":
        index_cmd = "minimap2 -d " + refname + " " + reference + " > /dev/null 2>&1"
    else:
        raise ValueError("Unrecognised mapper requested")

    exitcode = 0
    if not dry_run:
        exitcode = os.system(index_cmd)

    return index_cmd, exitcode


def map_reads(reads, reference, index, out_file, read_length, mapper="bwa", threads=1, dry_run=False):

    if mapper == "smalt":
        align_cmd = "smalt map -x -n " + str(threads) + " " + index + " " + reads + " 1> " + out_file + " 2> align.stderr"
    elif mapper == "bwa":
        k = min_seed_len_default(read_length)
        align_cmd = "bwa mem -k " + str(k) + " -t " + str(threads) + " " + reference + " " + reads + " 1> " + out_file + " 2> align.stderr"
    elif mapper == "minimap2":
        align_cmd = "minimap2 -c -o " + out_file + " -N 1 -ax sr " + reference + " " + reads
    elif mapper == "minimap2_long":
        align_cmd = "minimap2 -c -o " + out_file + " -N 1 -ax map-ont " + reference + " " + reads
    else:
        raise ValueError("Unrecognised mapper requested")

    exitcode = 0
    if not dry_run:
        exitcode = os.system(align_cmd)

    return align_cmd, exitcode


def min_seed_len_default(read_length):
    return 13 if read_length < 100 else 19


def smalt_k_default(read_length):
    return 13 if read_length < 100 else 20


def smalt_s_default(read_length):
    return 4 if read_length < 70 else 13 if read_length > 100 else 6


def sam2bam(sam_file_in, bam_file_out, threads=1):
    pysam.sort("-@", str(threads), "-O", "bam", "-T", sam_file_in + ".tmp", "-o", bam_file_out, sam_file_in)
    pysam.index("-@", str(threads), bam_file_out)
    stats = pysam.stats("-@", str(threads), bam_file_out)
    with open(bam_file_out + ".bamcheck", "w") as statsfile:
        statsfile.write(stats)
