"""
Combining insert site plot files
"""
import csv
import gzip
import os
import subprocess

from Bio import bgzf


class PlotFiles:

    def __init__(self, id, plot_files, working_dir):
        self.id = id
        self.plot_files = plot_files
        self.working_dir = working_dir

    @staticmethod
    def check_plot_files_lengths(plot_files, working_dir):
        lengths=[]
        for f in plot_files:
            plot_file = working_dir + os.sep + f if working_dir else f
            cat = "zcat" if plot_file[-3:] == ".gz" else "cat"
            wc = subprocess.check_output(["bash", "-c", cat + " " + plot_file + " | wc | awk '{print $1}'"])
            l = str(wc, 'UTF-8').strip()
            lengths.append(int(l))

        length = lengths[0]
        for x in range(1, len(lengths)):
            if lengths[x] != length:
                length = -1

        return length + 1, lengths;

    @staticmethod
    def create_file_handles(plot_files, working_dir=""):
        file_handles = []
        for f in plot_files:
            if f.endswith('.gz'):
                input_opener = gzip.open
                input_mode = "rt"
            else:
                input_opener = open
                input_mode = "r"

            file_handles.append(input_opener(working_dir + os.sep + f if working_dir else f, input_mode))
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
        for l in lines:
            if not l:
                return ""

            cols = l.strip().split(" ")
            forward_count = int(cols[0])
            reverse_count = int(cols[1])
            totals[0] += forward_count
            totals[1] += reverse_count

        return " ".join([str(x) for x in totals]), totals[0] or totals[1]

    @staticmethod
    def close_file_handles(file_handles):
        for f in file_handles:
            f.close()


def stats_header():
    return ["ID", "Sequence Length", "Unique Insertion Sites", "Seq Len/UIS"]


def get_plot_ids(plot_list_file):
    id_order = []
    with open(plot_list_file, 'r') as plf_handle:
        tsv_file = csv.reader(plf_handle, delimiter="\t")
        for row in tsv_file:
            id_order.append(row[0])
    return id_order


def prepare_and_create_tabix_for_combined_plots(tabix_plot, combined_dir):

    tabix_plot_name = combined_dir + os.sep + "tabix.insert_site_plot.gz"
    sorted_tabix_plot_name = combined_dir + os.sep + "tabix_sorted.insert_site_plot.gz"

    with bgzf.BgzfWriter(tabix_plot_name, "wb") as tabix_plot_fh:
        for line in tabix_plot:
            tabix_plot_fh.write(line + "\n")

    os.system("zcat " + tabix_plot_name + " | sort -k1,1 -k2,2n | bgzip > " + sorted_tabix_plot_name +
              " && tabix -b 2 -e 2 " + sorted_tabix_plot_name)
    os.remove(tabix_plot_name)


def combine(plot_file_list, combined_dir="combined"):

    with open(plot_file_list, 'r') as plf_handle:

        stats_file = os.path.splitext(os.path.basename(plot_file_list))[0] + ".stats"
        with open(stats_file, 'w') as stats_file_handle:

            os.makedirs(combined_dir, exist_ok=True)

            stats_writer = csv.writer(stats_file_handle, lineterminator="\n")
            stats_writer.writerow(stats_header())

            plf_reader = csv.reader(plf_handle, delimiter="\t")

            tabix_plot = []
            for line in plf_reader:
                pf = PlotFiles(line[0], line[1:], os.path.dirname(plot_file_list))
                tabix_data, length, nb_uis = pf.combine_plot_files(combined_dir)
                tabix_plot.extend(tabix_data)

                sl_per_uis = str(length / nb_uis) if nb_uis > 0 else "NaN"
                stats_writer.writerow([line[0], length, nb_uis, sl_per_uis])

            if len(tabix_plot) > 0:
                prepare_and_create_tabix_for_combined_plots(tabix_plot, combined_dir)
