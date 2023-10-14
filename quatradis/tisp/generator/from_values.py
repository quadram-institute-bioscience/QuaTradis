import numpy
import pyximport

pyximport.install()

from quatradis.tisp.write import write_plot_file, write_gzipped_plot_file


class PlotFromValuesGenerator:
    """
    Takes in arrays for forward and reverse integers and creates a new file
    """

    def __init__(self, forward, reverse, filename, gzipped=False):
        self.forward = forward
        self.reverse = reverse
        self.filename = filename
        self.gzipped = gzipped

        self.forward_length = len(self.forward)
        self.reverse_length = len(self.reverse)

    def construct_file(self):

        p_len = self.plot_length()
        if len(self.forward) == 0:
            self.forward = numpy.zeros(p_len, dtype=float)

        if len(self.reverse) == 0:
            self.reverse = numpy.zeros(p_len, dtype=float)

        if self.gzipped:
            write_gzipped_plot_file(self.filename, self.forward, self.reverse, p_len)
        else:
            write_plot_file(self.filename, self.forward, self.reverse, p_len)

        return self

    def plot_length(self):
        total_length = self.forward_length
        if self.reverse_length > total_length:
            total_length = self.reverse_length

        return total_length
