import os
import numpy
from quatradis.util.file_handle_helpers import ensure_output_dir_exists
from quatradis.tisp.parser import PlotParser
from quatradis.tisp.generator.from_values import PlotFromValuesGenerator


class NormalisePlots:
    '''Given a set of plot files, find the one with the highest number of reads and normalise all the rest into new files'''

    def __init__(self, plotfiles, minimum_proportion_insertions, output_dir=".", output_filename="normalised.plot_gz", verbose=False):
        self.plotfiles = plotfiles
        self.output_dir = output_dir
        self.output_filename = output_filename
        self.verbose = verbose
        self.minimum_proportion_insertions = minimum_proportion_insertions
        self.plot_objs = self.read_plots()

    def create_normalised_files(self):
        plot_objs, max_plot_reads = self.normalise()
        output_files = []
        for p in self.plotfiles:
            plotname = os.path.basename(p).split('.')[0]
            out_dir = os.path.join(self.output_dir, plotname)
            output_filename = os.path.join(out_dir, self.output_filename)
            ensure_output_dir_exists(out_dir)
            pg = PlotFromValuesGenerator(plot_objs[p].forward, plot_objs[p].reverse, output_filename)
            pg.construct_file()
            output_files.append(output_filename)
        # print("output_files: ",output_files )
        return output_files, max_plot_reads

    def decreased_insertion_reporting(self, max_plot_reads=0):
        if max_plot_reads == 0:
            max_plot_reads = self.max_reads(self.plot_objs)

        max_plot_insertions = self.max_insertions(self.plot_objs)

        normalising = True
        for p in self.plot_objs:
            t = self.plot_objs[p].total_reads
            if t / max_plot_reads < self.minimum_proportion_insertions:
                print("Number of reads in " + p + " is " + str(t) + " compared to a maximum of " + str(
                    max_plot_reads) + " so we cant call decreased insertions accurately")
                normalising = False

            t = self.plot_objs[p].total_insertions
            if t / max_plot_insertions < self.minimum_proportion_insertions:
                print("Number of insertions in " + p + " is " + str(t) + " compared to a maximum of " + str(
                    max_plot_insertions) + " so we cant call decreased insertions accurately")
                normalising = False
        return normalising

    def normalise(self):
        max_plot_reads = self.max_reads(self.plot_objs)
        print("max_plot_reads",max_plot_reads)

        for p in self.plotfiles:
            current_plot_reads = self.plot_objs[p].total_reads
            scaling_factor = max_plot_reads / current_plot_reads
            print("Plot file",p)
            print("scaling_factor",scaling_factor)
            print("****************")
            if self.verbose:
                print("\t".join(("Normalise", p, str(current_plot_reads), str(max_plot_reads), str(scaling_factor))))

            self.plot_objs[p].forward = numpy.multiply(self.plot_objs[p].forward, scaling_factor, dtype=float,
                                                       casting='unsafe')
            self.plot_objs[p].reverse = numpy.multiply(self.plot_objs[p].reverse, scaling_factor, dtype=float,
                                                       casting='unsafe')

        return self.plot_objs, max_plot_reads

    def read_plots(self):
        plot_objs = {}
        for p in self.plotfiles:
            pp = PlotParser(p)
            plot_objs[p] = pp
        return plot_objs

    def plot_insertions(self, plot_objs):
        ins = [plot_objs[p].total_insertions for p in plot_objs]
        return ins

    def plot_total_reads(self, plot_objs):
        reads = [plot_objs[p].total_reads for p in plot_objs]
        # print("plot_total_reads",reads)
        return reads

    def max_reads(self, plot_objs):
        return max(self.plot_total_reads(plot_objs))

    def max_insertions(self, plot_objs):
        return max(self.plot_insertions(plot_objs))
