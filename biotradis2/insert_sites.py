#!/usr/bin/env python3
# coding: utf-8

"""
Plotting insert sites
"""
import collections

import pysam
from Bio import bgzf


class SequencePlotter:
    def __init__(self, output_prefix, sequence_name, sequence_length):
        self.sequence_name = sequence_name
        self.cleaned_sequence_name = self.clean(sequence_name)
        self.output_file_handle = bgzf.BgzfWriter(output_prefix + "." + self.cleaned_sequence_name + ".insert_site_plot.gz", "wb")
        self.read_starts = collections.defaultdict(collections.Counter)
        self.sequence_length = sequence_length
        self.current_position = 0

    # Work out if padding is needed and return it as a formatted string.  Also update the current position.
    def create_padding_string(self, new_position):
        padding_string = ""
        for i in range(self.current_position + 1, new_position):
            padding_string += "0 0\n"
        self.current_position = new_position
        return padding_string

    # Returns the number of reads found at the specified position and direction
    def number_of_reads(self, read_coord, direction):
        if direction in self.read_starts[read_coord]:
            return self.read_starts[read_coord][direction]
        return 0

    # Prints out the current coordinate counts for both strands.  Automatically applies any necessary padding to file
    # prior to this read coord in necessary
    def print_coord(self, read_coord):
        padding_string = self.create_padding_string(read_coord)
        forward_reads = self.number_of_reads(read_coord, 1)
        reverse_reads = self.number_of_reads(read_coord, -1)
        self.output_file_handle.write(padding_string + forward_reads + " " + reverse_reads + "\n")

    def process_read_starts(self):
        for read_coord in sorted(self.read_starts.keys()):
            self.print_coord(read_coord)

        # Finish up this sequence with padding then close
        padding_string = self.create_padding_string(self.sequence_length)
        self.output_file_handle.write(padding_string)
        self.output_file_handle.close()



def count_read_starts(mapped_reads, mapping_score, plot_out):

    read_starts = {}

    # Samtools command to return only those alignments that are mapped and exceed the mapping score cutoff
    # samtools_command = "samtools view -F 4 -q " + str(mapping_score) + " " + mapped_reads
    samfile = pysam.AlignmentFile(mapped_reads, "rb")
    header = samfile.header
    for read in samfile.fetch():
        if not read.is_unmapped and read.mapping_quality >= mapping_score:

            sequence_name = read.get_reference_name()
            if sequence_name not in read_starts:
                read_starts[sequence_name] = SequencePlotter(plot_out, sequence_name, header.get_reference_length(sequence_name))

            sequence_plotter = read_starts[sequence_name]

            # TODO Might need to implement the cigar parser here to get proper start positions
            if read.is_reverse:
                sequence_plotter.read_starts[read.query_alignment_start][-1] += 1
            else:
                sequence_plotter.read_starts[read.query_alignment_start][1] += 1

    return read_starts


def plot(mapped_reads, plot_out, cutoff_score=30):
    """
    Generate insertion plots from a mapped fastq file
    :param mapped_reads: mapped and sorted BAM file
    :param plot_out: base name to assign to the resulting insertion site plot. Default = tradis.plot
    :param cutoff_score: cutoff value for mapping score. Default = 30
    :return:
    """

    seq_plotters = count_read_starts(mapped_reads, cutoff_score, plot_out)
    for sequence_name in sorted(seq_plotters.keys()):
        seq_plotters[sequence_name].process_read_starts()

