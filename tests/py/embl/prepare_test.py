import unittest
import os
import filecmp

from quatradis.embl.prepare import PrepareEMBLFile
from tests.py.embl import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'prepareinputfiles')
data_expand_genes_dir = os.path.join(DATA_DIR, 'expandgenes')


class TestPrepareEMBLFile(unittest.TestCase):

	def test_small_valid_file(self):
		p = PrepareEMBLFile(os.path.join(data_dir, 'valid'), 1, 4, 2, 100, None)
		
		self.assertTrue(p.create_file())
		self.assertTrue(os.path.exists(p.embl_filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_embl_filename'), p.embl_filename))

		os.remove(p.embl_filename)	
		
	def test_reuse_annotation(self):
		p = PrepareEMBLFile(os.path.join(data_dir, 'valid'), 1, 4, 2, 15, os.path.join(data_expand_genes_dir, 'one_gene'))
		self.assertTrue(p.create_file())
		
		self.assertTrue(filecmp.cmp(os.path.join(data_expand_genes_dir, 'expected_one_gene'), p.embl_filename))
		os.remove(p.embl_filename)
