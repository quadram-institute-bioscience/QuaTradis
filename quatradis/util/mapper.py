"""
Mapping prepared reads to index
"""
import gzip
import os

import pysam
from Bio import SeqIO

from quatradis.util.file_handle_helpers import reader_opener

def index_and_align(reads, reference, index, out_file, mapper="bwa", threads=1):
    read_length = calc_read_length(reads)
    index_cmd = index_reference(reference, index, read_length, mapper)
    align_cmd = map_reads(reads, reference, index, out_file, read_length, mapper, threads)
    return index_cmd, align_cmd

def extract_sequence_names(reference, out_file):
    with open(out_file, "w") as f:
        input_opener, input_mode = reader_opener(reference)
        read_length = 0
        with input_opener(reference, input_mode) as inputHandle:
            for seq_record in SeqIO.parse(inputHandle, "fasta"):
                f.write(str(seq_record.id) + "\n")

def calc_read_length(reads_file):

    input_opener, input_mode = reader_opener(reads_file)
    read_length = 0
    with input_opener(reads_file, input_mode) as inputHandle:
        records = SeqIO.parse(inputHandle, "fastq")
        for rec in records:
            read_length = len(rec.seq)
            break

    inputHandle.close()

    return read_length


def index_reference(reference, refname, read_length, mapper="bwa", dry_run=False):
    ignore_output = " > /dev/null 2>&1"
    if mapper == "smalt":
        k = smalt_k_default(read_length)
        s = smalt_s_default(read_length)
        index_cmd = "smalt index -k " + str(k) + " -s " + str(s) + " " + refname + " " + reference + ignore_output
    elif mapper == "bwa":
        index_cmd = "bwa index " + reference + ignore_output
    elif mapper == "minimap2" or mapper == "minimap2_long":
        index_cmd = "minimap2 -d " + refname + " " + reference + ignore_output
    else:
        raise ValueError("Unrecognised mapper requested")

    exitcode = 0
    if not dry_run:
        exitcode = os.system(index_cmd)

    return index_cmd, exitcode


def map_reads(reads, reference, index, out_file, read_length, mapper="bwa", threads=1, dry_run=False):

    if mapper == "smalt":
        align_cmd = "smalt map -x -n " + str(threads) + " " + index + " " + reads + " 1> " + out_file + " 2> " + out_file + ".stderr"
    elif mapper == "bwa":
        k = min_seed_len_default(read_length)
        align_cmd = "bwa mem -k " + str(k) + " -t " + str(threads) + " " + reference + " " + reads + " 1> " + out_file + " 2> " + out_file + ".stderr"
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
    if read_length < 70:
        return 4
    elif read_length > 100:
        return 13
    else:
        return 6


def sam2bam(sam_file_in, bam_file_out, threads=1):
    pysam.sort("-@", str(threads), "-O", "bam", "-T", sam_file_in + ".tmp", "-o", bam_file_out, sam_file_in)
    pysam.index("-@", str(threads), bam_file_out)
    stats = pysam.stats("-@", str(threads), bam_file_out)
    with open(bam_file_out + ".bamcheck", "w") as statsfile:
        statsfile.write(stats)
