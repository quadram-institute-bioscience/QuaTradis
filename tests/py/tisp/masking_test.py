import unittest
import os
import filecmp

from quatradis.tisp.masking import PlotMasking
from tests.py.tisp import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'masking')


class TestPlotMasking(unittest.TestCase):
	
	def test_plot_masking_small(self):
		case = os.path.join(data_dir, 'sample1.plot')
		control = os.path.join(data_dir, 'sample2.plot')
		mask = os.path.join(data_dir, 'mask.plot')
		pm = PlotMasking([case, control], mask, True)
		self.assertEqual(2, len(pm.output_plot_files))
		
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_sample1.plot'), pm.output_plot_files[case]))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_sample2.plot'), pm.output_plot_files[control]))

	def test_no_masking(self):
		case = os.path.join(data_dir, 'sample1.plot')
		control = os.path.join(data_dir, 'sample2.plot')
		mask = os.path.join(data_dir, 'mask.plot')
		pm = PlotMasking([case, control], mask, False)
		self.assertEqual(2, len(pm.output_plot_files))
		
		self.assertTrue(filecmp.cmp(case, pm.output_plot_files[case]))
		self.assertTrue(filecmp.cmp(control, pm.output_plot_files[control]))
