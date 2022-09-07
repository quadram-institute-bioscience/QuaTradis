from tempfile import mkstemp

from quatradis.tisp.parser import PlotParser

from quatradis.embl.window import WindowGenerator
from quatradis.embl.generator import EMBLGenerator
from quatradis.embl.expand_genes import EMBLExpandGenes


class PrepareEMBLFile:
    '''Take in the input files, parse them to create new files.'''

    def __init__(self, plotfile, minimum_threshold, window_size, window_interval, prime_feature_size, emblfile):
        self.plotfile = plotfile
        self.minimum_threshold = minimum_threshold
        self.window_size = window_size
        self.window_interval = window_interval
        self.prime_feature_size = prime_feature_size
        self.emblfile = emblfile

        self.forward_plot_filename = ""
        self.reverse_plot_filename = ""
        self.combined_plot_filename = ""
        self.embl_filename = ""

    def genome_length(self):
        return self.plot_parser_obj.genome_length

    def plot_parser(self):
        return PlotParser(self.plotfile, self.minimum_threshold)

    def windows(self):
        w = WindowGenerator(self.plot_parser_obj.genome_length, self.window_size, self.window_interval)
        return w.create_windows()

    def create_embl_file(self, embl_filename):
        e = EMBLGenerator(self.windows(), self.plot_parser_obj.genome_length)

        if not embl_filename:
            fd, embl_filename = mkstemp()
        e.construct_file(embl_filename)
        return embl_filename

    def embl_file_expand_genes(self, embl_filename):
        if not embl_filename:
            fd, embl_filename = mkstemp()

        eg = EMBLExpandGenes(self.emblfile, self.prime_feature_size)
        eg.construct_file(embl_filename)
        return embl_filename

    def create_file(self, embl_filename=None):
        self.plot_parser_obj = self.plot_parser()

        # use the annotation and add 3/5 prime blocks to each gene
        if self.emblfile:
            self.embl_filename = self.embl_file_expand_genes(embl_filename)
        else:
            # Use sliding windows only
            self.embl_filename = self.create_embl_file(embl_filename)
        return self.embl_filename
