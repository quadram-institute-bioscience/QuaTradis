import numpy
import os
import pandas
import subprocess

from Bio import bgzf

import quatradis.util.file_handle_helpers as fhh


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass


class PlotFiles:

    def __init__(self, id, plot_files, working_dir):
        self.id = id
        self.plot_files = plot_files
        self.working_dir = working_dir

    @staticmethod
    def check_plot_files_lengths(plot_files, working_dir):
        lengths = []
        for f in plot_files:
            lengths.append(PlotParser.get_plot_file_length(f, working_dir))

        length = lengths[0]
        for x in range(1, len(lengths)):
            if lengths[x] != length:
                length = -1

        return length + 1, lengths

    @staticmethod
    def create_file_handles(plot_files, working_dir=""):
        file_handles = []
        for f in plot_files:
            file_handles.append(PlotParser.create_file_handle(f, working_dir))
        return file_handles

    def combine_plot_files(self, combined_dir):

        length, lengths = PlotFiles.check_plot_files_lengths(self.plot_files, self.working_dir)
        if length == -1:
            raise ValueError("Inconsistent file lengths for \"" + self.id + "\" : " + lengths)

        file_handles = self.create_file_handles(self.plot_files, self.working_dir)

        tabix_plot = []
        full_plot = []
        nb_uis = 0
        for i in range(length):
            combined_line, uis = self.combine_next_lines(file_handles, i + 1)
            if uis:
                nb_uis += 1

            plot_values_tabix = combined_line.replace(' ', '\t')

            tabix_line = "\t".join([self.id, str(i), plot_values_tabix])
            tabix_plot.append(tabix_line)
            full_plot.append(combined_line)

        self.close_file_handles(file_handles)

        comb_plot_name = combined_dir + os.sep + self.id + ".insert_site_plot.gz"
        with bgzf.BgzfWriter(comb_plot_name, "wb") as cplot:
            for line in full_plot:
                cplot.write(line + "\n")

        return tabix_plot, length, nb_uis

    def combine_next_lines(self, file_handles, line_number):
        current_lines = []
        for current_handle in file_handles:
            line = current_handle.readline()
            if line:
                current_lines.append(line)
            else:
                self.close_file_handles(file_handles)
                raise ValueError("Unexpected blank line in \"" + self.id + "\" line " + (line_number))

        return PlotFiles.combine_lines(current_lines)

    @staticmethod
    def combine_lines(lines):
        """
        Combine plot file lines at same base position from multiple samples
        :param lines: List of plot file lines
        :return: Summed plot file line, and whether or not this combined base position is a unique insertion site
        """
        totals = [0, 0]
        for line in lines:
            if not line:
                return ""

            cols = line.strip().split(" ")
            forward_count = int(cols[0])
            reverse_count = int(cols[1])
            totals[0] += forward_count
            totals[1] += reverse_count

        return " ".join([str(x) for x in totals]), totals[0] or totals[1]

    @staticmethod
    def close_file_handles(file_handles):
        for f in file_handles:
            f.close()


class PlotParser:

    def __init__(self, filename, minimum_threshold=0):
        self.filename = filename
        self.minimum_threshold = minimum_threshold
        self.forward = []
        self.reverse = []
        self.combined = []
        self.genome_length = 0

        self.read()

        self.total_reads = sum(self.combined)
        self.total_insertions = sum([1 for a in self.combined if a > 0])

    @staticmethod
    def create_file_handle(plot_file, working_dir=""):
        input_opener, input_mode = fhh.reader_opener(plot_file)
        return input_opener(os.path.join(working_dir, plot_file) if working_dir else plot_file, input_mode)

    @staticmethod
    def get_plot_file_length(plot_file, working_dir):
        plot_file = working_dir + os.sep + plot_file if working_dir else plot_file
        cat = "zcat" if plot_file[-3:] == ".gz" else "cat"
        wc = subprocess.check_output(["bash", "-c", cat + " " + plot_file + " | wc | awk '{print $1}'"])
        file_length = str(wc, 'UTF-8').strip()
        return int(file_length)

    def read(self):
        handle = PlotParser.create_file_handle(self.filename)
        insert_site_array = pandas.read_csv(handle, sep='\s+', dtype=float, engine='c',
                                            header=None).values

        self.genome_length = len(insert_site_array)

        self.forward = insert_site_array[:, 0]
        self.reverse = insert_site_array[:, 1]

        if self.minimum_threshold != 0:
            self.forward = self.filter_column(self.forward, self.genome_length)
            self.reverse = self.filter_column(self.reverse, self.genome_length)

        self.combined = [self.forward[i] + self.reverse[i] for i in range(0, self.genome_length)]

        return self

    def filter_column(self, ins_array, genome_length):
        abs_ins_values = numpy.absolute(ins_array)
        return [ins_array[i] if abs_ins_values[i] >= self.minimum_threshold else 0 for i in range(0, genome_length)]

    def __str__(self):
        return "\t".join((self.filename, str(len(self.combined)), str(self.total_reads), str(self.total_insertions),
                          str(self.total_reads / self.total_insertions)))


class ScoreParser:
    def __init__(self, filename, minimum_threshold=0):
        self.filename = filename
        self.minimum_threshold = minimum_threshold
        self.pvals = []
        self.qvals = []
        self.genome_length = 0

        self.split_lines()

    def split_lines(self):
        insert_site_array = pandas.read_csv(self.filename, delim_whitespace=True, dtype=float, engine='c',
                                            header=None).values

        self.genome_length = len(insert_site_array)

        self.pvals = insert_site_array[:, 0]
        self.qvals = insert_site_array[:, 1]

        if self.minimum_threshold != 0:
            self.pvals = self.filter_column(self.pvals, self.genome_length)
            self.qvals = self.filter_column(self.qvals, self.genome_length)

        return self

    def filter_column(self, ins_array, genome_length):
        abs_ins_values = numpy.absolute(ins_array)
        return [ins_array[i] if abs_ins_values[i] >= self.minimum_threshold else 0 for i in range(0, genome_length)]
