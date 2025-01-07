"""
Creating insert site plot files (originally designed for artemis)
"""
import collections
import os
import subprocess
import sys

import pysam
from Bio import bgzf
from pysam.libcalignedsegment import CIGAR_OPS

from quatradis.util.file_handle_helpers import input_alignment_mode


class PlotFromAlignmentsGenerator:
    def __init__(self, output_prefix, sequence_name, sequence_length):
        self.sequence_name = sequence_name
        self.cleaned_sequence_name = self.clean_sequence_name(sequence_name)
        if os.path.dirname(output_prefix):
            os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        self.output_file_handle = bgzf.BgzfWriter(output_prefix + "." + self.cleaned_sequence_name + ".insert_site_plot.gz", "wb")
        self.read_starts = collections.defaultdict(collections.Counter)
        self.sequence_length = sequence_length
        self.current_position = 0

    def clean_sequence_name(self, sequence_name):
        # TODO Make sure sequence name doesn't contain symbols or whitespace etc.
        return sequence_name

    # Work out if padding is needed and return it as a formatted string.  Also update the current position.
    def create_padding_string(self, new_position):
        padding_string = ""
        for _ in range(self.current_position + 1, new_position):
            padding_string += "0 0\n"
        self.current_position = new_position
        return padding_string

    # Returns the number of reads found at the specified position and direction
    def number_of_reads(self, read_coord, direction):
        if direction in self.read_starts[read_coord]:
            return self.read_starts[read_coord][direction]
        return 0

    # Prints out the current coordinate counts for both strands.  Automatically applies any necessary padding to file
    # prior to this read coord if necessary
    def print_coord(self, read_coord):
        padding_string = self.create_padding_string(read_coord)
        forward_reads = self.number_of_reads(read_coord, 1)
        reverse_reads = self.number_of_reads(read_coord, -1)
        self.output_file_handle.write(padding_string + str(forward_reads) + " " + str(reverse_reads) + "\n")

    def process_read_starts(self):
        """
        Actually creates the output plot file by processing all read start positions in order.  Before writing the
        count at the current position it pads any skipped positions with '0 0'.

        There is a hack implemented here to avoid accidentally creating plot files that are too large.
        Any reads that start off the end of the reference, which can happen from alignments to circular sequences,
        are ignored.  We just count how often we are ignoring them and return that value along with the number of
        insertion sites found in this reference sequence.
        """
        nb_unique_insertion_sites = 0
        nb_skipped = 0
        for read_coord in sorted(self.read_starts.keys()):
            if read_coord > self.sequence_length:
                nb_skipped += 1
                continue
            self.print_coord(read_coord)
            nb_unique_insertion_sites += 1

        # Finish up this sequence with any required padding then close the file
        padding_string = self.create_padding_string(self.sequence_length+1)
        self.output_file_handle.write(padding_string)
        self.output_file_handle.close()

        return nb_unique_insertion_sites, nb_skipped


def op_is_aligned(op):
    return op == CIGAR_OPS.CMATCH or op == CIGAR_OPS.CDIFF or op == CIGAR_OPS.CEQUAL


def op_is_not_on_ref(op):
    return op == CIGAR_OPS.CSOFT_CLIP or op == CIGAR_OPS.CDEL or op == CIGAR_OPS.CREF_SKIP


def process_op(op, number, start, end, current_coordinate):
    if op_is_aligned(op):
        if start == 0:
            start = current_coordinate
        current_coordinate += number
        if end < current_coordinate:
            end = current_coordinate - 1
    elif op_is_not_on_ref(op):
        if start == 0:
            current_coordinate -= number
            start = current_coordinate
        current_coordinate += number
        if end < current_coordinate:
            end = current_coordinate - 1
    # Ignore inserts
    #elif op == CIGAR_OPS.CINS:
    #    continue

    return start, end, current_coordinate


def cigar_parser(cigar, coord):
    start = 0
    end = 0
    current_coordinate = coord
    for action, number in cigar:
        op = CIGAR_OPS(action)
        start, end, current_coordinate = process_op(op, number, start, end, current_coordinate)

    return start, end


def count_read_starts(mapped_reads, cutoff_score, plot_out_prefix):

    read_starts = {}
    nb_mapped = 0
    nb_mapped_and_good = 0

    # Samtools command to return only those alignments that are mapped and exceed the mapping score cutoff
    # samtools_command = "samtools view -F 4 -q " + str(mapping_score) + " " + mapped_reads
    alignments = pysam.AlignmentFile(mapped_reads, input_alignment_mode(mapped_reads))

    for ref in alignments.header.references:
        read_starts[ref] = PlotFromAlignmentsGenerator(plot_out_prefix, ref, alignments.header.get_reference_length(ref))

    for read in alignments.fetch():
        if not read.is_unmapped:
            nb_mapped += 1
            if read.mapping_quality >= cutoff_score:
                nb_mapped_and_good += 1
                sequence_plotter = read_starts[read.reference_name]
                start, end = cigar_parser(read.cigartuples, read.pos+1)
                if read.is_reverse:
                    sequence_plotter.read_starts[end][-1] += 1
                else:
                    sequence_plotter.read_starts[start][1] += 1

    return read_starts, nb_mapped, nb_mapped_and_good


def create_plot_files(mapped_reads, plot_out_prefix="tradis.plot", cutoff_score=30):

    uis_map = {}
    is_plotters, nb_mapped, nb_mapped_and_good = count_read_starts(mapped_reads, cutoff_score, plot_out_prefix)
    for sequence_name in sorted(is_plotters.keys()):
        nb_unique_insertion_sites, nb_skipped = is_plotters[sequence_name].process_read_starts()
        uis_map[sequence_name] = nb_unique_insertion_sites
        if nb_skipped > 0:
            print("Skipped", nb_skipped, "insertion sites which were found off the end of the reference sequence:", sequence_name)

    return nb_mapped, nb_mapped_and_good, uis_map


def get_number_reads(seq_file_in):

    # If running on mac take from stdin
    joiner = " < " if sys.platform == "darwin" else " "

    # Make sure we can handle gzipped or non-gzipped fastq input
    if seq_file_in.endswith('.gz'):
        cat_cmd = "zcat"
    else:
        cat_cmd = "cat"
    lines = int(subprocess.check_output(["bash", "-c", cat_cmd + joiner + seq_file_in + " | wc -l"]).strip())
    return lines / 4


def generate_stats(mapped_reads, out_file, fastq, nb_reads, nb_tagged_reads, nb_mapped, uis_map):

    with open(out_file, "w") as output_handle:

        samfile = pysam.AlignmentFile(mapped_reads, "rb")
        sequences = samfile.header.references

        header = create_stats_file_header(sequences)
        output_handle.write(",".join(header) + "\n")

        output_filename = os.path.split(fastq)[1]

        if nb_reads == 0:
            nb_reads = get_number_reads(fastq)

        if nb_tagged_reads == 0:
            nb_tagged_reads = nb_reads

        matching_pc = (nb_tagged_reads / nb_reads) * 100.0
        mapped_pc = (nb_mapped / nb_tagged_reads) * 100.0
        row = [output_filename, str(nb_reads), str(nb_tagged_reads), str(matching_pc), str(nb_mapped), str(mapped_pc)]

        total_unique_insertion_sites = 0
        total_seq_length = 0
        for sn in sorted(sequences):
            unique_insertion_sites = uis_map[sn]
            total_unique_insertion_sites += unique_insertion_sites
            sequence_length = samfile.header.get_reference_length(sn)
            total_seq_length += sequence_length
            uis_per_seqlen = str(sequence_length / unique_insertion_sites) if unique_insertion_sites > 0 else "NaN"
            row.extend([str(unique_insertion_sites), uis_per_seqlen])

        total_uis_per_total_seqlen = str(total_seq_length / total_unique_insertion_sites) if total_unique_insertion_sites > 0 else "NaN"
        row.extend([str(total_unique_insertion_sites), str(total_uis_per_total_seqlen)])

        output_handle.write(",".join(row))


def create_stats_file_header(sequence_names):

    fields = [
        "File",
        "Total Reads",
        "Reads Matched",
        "% Matched",
        "Reads Mapped",
        "% Mapped"]

    for sn in sequence_names:
        fields.append("Unique Insertion Sites : " + sn)
        fields.append("Seq Len/UIS : " + sn)

    fields.append("Total Unique Insertion Sites")
    fields.append("Total Seq Len/Total UIS")

    return fields


def plot(mapped_reads, fastq="", plot_out_prefix="tradis.plot", cutoff_score=30, nb_reads=0, nb_tagged_reads=0):
    """
    Generate insertion plots from a mapped fastq file and calculate stats
    :param mapped_reads: mapped and sorted BAM file
    :param fastq: the original fastq file used as input for the mapper
    :param plot_out_prefix: base name to assign to the resulting insertion site plot. Default = tradis.plot
    :param nb_reads: the original number of reads in the fastq file.  If set to 0, we will run through the fastq to count them.  Default=0
    :param nb_tagged_reads: the number of reads are tags are identified and clipped.  If set to 0 we assume it is the same as nb_reads.  Default 0.
    :param cutoff_score: cutoff value for mapping score. Default = 30
    :return:
    """

    # Create the plot files
    nb_mapped, nb_mapped_and_good, uis_map = create_plot_files(mapped_reads, plot_out_prefix, cutoff_score)

    if fastq:
        # Calculate and output the stats in a csv
        generate_stats(mapped_reads, plot_out_prefix + ".stats", fastq, nb_reads, nb_tagged_reads, nb_mapped, uis_map)
