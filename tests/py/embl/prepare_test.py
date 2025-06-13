import unittest
import os
import filecmp

from quatradis.embl.prepare import PrepareEMBLFile
from quatradis.util.config_defaults import DYNAMIC_WINDOW_PARAMS
from tests.py.embl import DATA_DIR
 
data_dir = os.path.join(DATA_DIR, 'prepareinputfiles')
data_expand_genes_dir = os.path.join(DATA_DIR, 'expandgenes')



class TestPrepareEMBLFile(unittest.TestCase):

	def test_small_valid_file(self):
		dynamic_window_params= {key: None for key, value in DYNAMIC_WINDOW_PARAMS.items()}
		dynamic_window_params.pop("dynamic_window")
		plotfiles= [os.path.join(data_dir, 'valid'),os.path.join(data_dir, 'valid_ctrl')]
		p = PrepareEMBLFile(plotfiles, 1, 4, 2, 15, None,False,False,dynamic_window_params)
		# p = PrepareEMBLFile(os.path.join(data_dir, 'valid'), 1, 4, 2, 100, None)
		self.assertTrue(p.create_file())
		self.assertTrue(os.path.exists(p.embl_filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_embl_filename'), p.embl_filename))
		os.remove(p.embl_filename)	
	
	#reuse annotation without dynamic window
	def test_reuse_annotation(self):
		dynamic_window_params= {key: None for key, value in DYNAMIC_WINDOW_PARAMS.items()}
		dynamic_window_params.pop("dynamic_window")
		plotfiles= [os.path.join(data_dir, 'valid'),os.path.join(data_dir, 'valid_ctrl')]
		p = PrepareEMBLFile(plotfiles, 1, 4, 2, 15, os.path.join(data_expand_genes_dir, 'one_gene'),False,False,dynamic_window_params)
		self.assertTrue(p.create_file())
		
		self.assertTrue(filecmp.cmp(os.path.join(data_expand_genes_dir, 'expected_one_gene'), p.embl_filename))
		os.remove(p.embl_filename)
    
	#reuse annotation with dynamic window
	def test_prepare_embl_dynamic_window(self):
		dynamic_window_params= {key: value[2] for key, value in DYNAMIC_WINDOW_PARAMS.items()}
		dynamic_window_params.pop("dynamic_window")
		plotfiles= [os.path.join(data_dir, 'valid'),os.path.join(data_dir, 'valid_ctrl')]
		p = PrepareEMBLFile(plotfiles, 1, 4, 2, 15, os.path.join(data_expand_genes_dir, 'one_gene'),True,False,dynamic_window_params)
		self.assertTrue(p.create_file())

