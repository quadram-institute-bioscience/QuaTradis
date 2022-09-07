import filecmp
import os
import unittest

from quatradis.tisp.generator.from_values import PlotFromValuesGenerator
from quatradis.tisp.parser import PlotParser

data_dir = os.path.join('data', 'tisp', 'generator', 'from_values')


class TestPlotGenerator(unittest.TestCase):

    def test_forward_only(self):
        filename = os.path.join(data_dir, 'test_plotgen')
        p = PlotFromValuesGenerator([0, 0, 0, 0, 0, 1, 1, 3], [], filename)
        self.assertTrue(p.construct_file())

        self.assertTrue(os.path.exists(filename))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_plotgen'), filename))
        os.remove(filename)

    def test_generate_and_read(self):
        filename = os.path.join(data_dir, 'plot.test')
        p = PlotFromValuesGenerator([0, 0, 0, 0, 0, 1, 1, 3], [9, 9, 0, 9, 9, 1, 1, 3], filename)
        p.construct_file()

        d = PlotParser(filename, 0)
        self.assertTrue(self.check_arrays_equal(d.forward, [0, 0, 0, 0, 0, 1, 1, 3]))
        self.assertTrue(self.check_arrays_equal(d.reverse, [9, 9, 0, 9, 9, 1, 1, 3]))
        self.assertTrue(self.check_arrays_equal(d.combined, [9, 9, 0, 9, 9, 2, 2, 6]))

        os.remove(filename)

    def check_arrays_equal(self, array1, array2):
        for i, val in enumerate(array1):
            if array1[i] != array2[i]:
                print(str(array1[i]) + " not equal to " + str(array2[i]))
                return False
        return True
