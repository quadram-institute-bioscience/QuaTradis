from tempfile import mkstemp

from quatradis.tisp.parser import PlotParser

from quatradis.embl.window import WindowGenerator
from quatradis.embl.generator import EMBLGenerator
from quatradis.embl.expand_genes import EMBLExpandGenes


class PrepareEMBLFile:
    '''Take in the input files, parse them to create new files.'''
    # Modification 4
    def __init__(self, plotfile, minimum_threshold, window_size, window_interval,prime_feature_size, emblfile,dynamic_window,disable_new_algorithm,kwargs):
        #First two would be condition files and last two would be control files.
        self.plotfiles = plotfile
        self.minimum_threshold = minimum_threshold
        self.window_size = window_size
        self.window_interval = window_interval
        self.prime_feature_size=prime_feature_size
        # Optional: Store kwargs as a dictionary for further reference
        self.dynamic_params = kwargs 
        self.emblfile = emblfile
        self.dynamic_window= dynamic_window
        self.disable_new_algorithm=disable_new_algorithm
        self.forward_plot_filename = ""
        self.reverse_plot_filename = ""
        self.combined_plot_filename = ""
        self.embl_filename = ""

    def genome_length(self):
        return self.plot_parser_objs["Control1"].genome_length

    def plot_parser(self):
        plot_parsers = {}
        num_files = len(self.plotfiles)

        if num_files % 2 != 0:
            raise ValueError("The number of plotfiles must be even.")

        half_index = num_files // 2

        for i in range(half_index):
            condition_file = self.plotfiles[i]
            control_file = self.plotfiles[i + half_index]
            plot_parsers[f"Condition{i+1}"] = PlotParser(condition_file,self.minimum_threshold)
            plot_parsers[f"Control{i+1}"] = PlotParser(control_file,self.minimum_threshold)

        return plot_parsers

    def windows(self):
        w = WindowGenerator(self.plot_parser_objs["Control1"].genome_length, self.window_size, self.window_interval)
        return w.create_windows()

    def create_embl_file(self, embl_filename):
        print("Genome Length (EMBL File Creation)",)
        e = EMBLGenerator(self.windows(), self.plot_parser_objs["Control1"].genome_length)

        if not embl_filename:
            fd, embl_filename = mkstemp()
        e.construct_file(embl_filename)
        return embl_filename

    def embl_file_expand_genes(self, embl_filename):
        if not embl_filename:
            fd, embl_filename = mkstemp()
        # Modification 5
        eg = EMBLExpandGenes(self.emblfile,self.prime_feature_size,self.dynamic_window,self.disable_new_algorithm,self.dynamic_params)
        eg.construct_file(embl_filename, self.plot_parser_objs)
        return embl_filename

    def create_file(self, embl_filename=None):
        self.plot_parser_objs = self.plot_parser()

        # use the annotation and add 3/5 prime blocks to each gene
        if self.emblfile:
            self.embl_filename = self.embl_file_expand_genes(embl_filename)
        else:
            # Use sliding windows only
            self.embl_filename = self.create_embl_file(embl_filename)
        return self.embl_filename
