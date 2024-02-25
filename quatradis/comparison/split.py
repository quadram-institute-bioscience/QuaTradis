import os

from quatradis.tisp.generator.from_values import PlotFromValuesGenerator
from quatradis.tisp.parser import PlotParser


class SplitPlotFile:
    """Used to take in a plot file and split it into forward only, reverse only and combined plot files."""

    def __init__(self, plotfile, minimum_threshold, gzipped=False, output_dir=None):
        self.plotfile = plotfile
        self.minimum_threshold = minimum_threshold
        self.gzipped = gzipped
        self.ext = ".gz" if self.gzipped else ""
        self.output_dir = "." if not output_dir else output_dir

    def _construct_file_path(self, type):
        return os.path.join(self.output_dir, type + ".plot" + self.ext)

    def get_forward_file_path(self):
        return self._construct_file_path("forward")

    def get_reverse_file_path(self):
        return self._construct_file_path("reverse")

    def get_combined_file_path(self):
        return self._construct_file_path("combined")

    def _create_split_plot_file(self, forward, reverse, filename):
        p = PlotFromValuesGenerator(forward, reverse, filename)
        p.construct_file()

    def create_all_files(self):
        plot_parser_obj = PlotParser(self.plotfile, self.minimum_threshold)

        os.makedirs(self.output_dir, exist_ok=True)

        self._create_split_plot_file(plot_parser_obj.forward, [], self.get_forward_file_path())
        self._create_split_plot_file([], plot_parser_obj.reverse, self.get_reverse_file_path())
        self._create_split_plot_file(plot_parser_obj.forward, plot_parser_obj.reverse, self.get_combined_file_path())
        return plot_parser_obj.genome_length


def split_plot(plot_file, output_dir, minimum_threshold=5, gzipped=True, verbose=False):
    p = SplitPlotFile(plot_file, minimum_threshold, gzipped=gzipped, output_dir=output_dir)
    length = p.create_all_files()
    if verbose:
        print("Split plot files.  Genome Length: " + str(length))
        print("Forward plot file: " + p.get_forward_file_path())
        print("Reverse plot file: " + p.get_reverse_file_path())
        print("Combined plot file: " + p.get_combined_file_path())
