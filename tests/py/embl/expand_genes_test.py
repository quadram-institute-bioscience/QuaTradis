import unittest
import os
import filecmp

from quatradis.embl.expand_genes import EMBLExpandGenes
from tests.py.embl import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'expandgenes')

class TestEMBLExpandGenes(unittest.TestCase):

	def test_embl_expand_genes(self):
		e = EMBLExpandGenes(os.path.join(data_dir, 'one_gene') , 15)
		
		output_file = os.path.join(data_dir, 'actual_one_gene')
		e.construct_file(output_file)
		self.assertTrue(os.path.exists(output_file))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_one_gene'), output_file))
		
		os.remove(output_file)	
		