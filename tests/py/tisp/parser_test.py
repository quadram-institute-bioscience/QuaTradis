import os
import numpy as np
import unittest


from quatradis.tisp.parser import PlotParser
from tests.py.tisp import DATA_DIR


class ErrorReadingFile (Exception): pass


class InvalidFileFormat (Exception): pass


data_dir = os.path.join(DATA_DIR, "parser")


class ParserPlotsTest(unittest.TestCase):

	def test_valid(self):
		tests = [
			('valid', 0, [0,0,0,0,0,1,1,3], [0,1,4,5,0,0,0,0], [0,1,4,5,0,1,1,3]),
			('valid.gz', 0, [0,0,0,0,0,1,1,3], [0,1,4,5,0,0,0,0], [0,1,4,5,0,1,1,3]),
			('valid', 3, [0,0,0,0,0,0,0,3], [0,0,4,5,0,0,0,0], [0,0,4,5,0,0,0,3])
		]
		for i, t in enumerate(tests):
			with self.subTest(i=i):
				d = PlotParser(os.path.join(data_dir, t[0]), t[1])
				assert np.array_equal(d.forward, t[2])
				assert np.array_equal(d.reverse, t[3])
				assert np.array_equal(d.combined, t[4])

	def test_large_file_zipped(self):
		d = PlotParser(os.path.join(data_dir,'Control2.out.CP009273.insert_site_plot.gz'),3)
		self.assertEqual(d.total_reads, 7684360)

