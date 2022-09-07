import os
from tempfile import mkstemp

from quatradis.tisp.parser import PlotParser
from quatradis.tisp.generator.from_values import PlotFromValuesGenerator


class PrepareInputFiles:
    '''Take in the input files, parse them to create new files.'''

    def __init__(self, plotfile, minimum_threshold, gzipped=False):
        self.plotfile = plotfile
        self.minimum_threshold = minimum_threshold
        self.gzipped = gzipped

        self.forward_plot_filename = ""
        self.reverse_plot_filename = ""
        self.combined_plot_filename = ""
        self.original_plot_filename = plotfile
        self.embl_filename = ""

    def genome_length(self):
        return self.plot_parser_obj.genome_length

    def plot_parser(self):
        return PlotParser(self.plotfile, self.minimum_threshold)

    def create_split_plot_file(self, forward, reverse, filename=None):
        if not filename:
            fd, filename = mkstemp()
        p = PlotFromValuesGenerator(forward, reverse, filename, self.gzipped)
        p.construct_file()
        return filename

    def create_all_files(self, output_dir=None):
        self.plot_parser_obj = self.plot_parser()

        ext = "plot"
        if self.gzipped:
            ext += ".gz"
        self.forward_plot_filename = self.create_split_plot_file(self.plot_parser_obj.forward, [], os.path.join(output_dir, "forward."+ext) if output_dir else None)
        self.reverse_plot_filename = self.create_split_plot_file([], self.plot_parser_obj.reverse, os.path.join(output_dir, "reverse."+ext) if output_dir else None)
        self.combined_plot_filename = self.create_split_plot_file(self.plot_parser_obj.forward, self.plot_parser_obj.reverse, os.path.join(output_dir, "combined."+ext) if output_dir else None)
        return self
