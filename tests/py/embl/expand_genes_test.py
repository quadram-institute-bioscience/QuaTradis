import unittest
import os
import filecmp

from quatradis.embl.expand_genes import EMBLExpandGenes
from quatradis.embl.prepare import PrepareEMBLFile
from quatradis.util.config_defaults import DYNAMIC_WINDOW_PARAMS
from tests.py.embl import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'expandgenes')

class TestEMBLExpandGenes(unittest.TestCase):
    
	#expand genes using dynamic window
	def test_embl_expand_genes_dw(self):
		dynamic_window_params= {key: value[2] for key, value in DYNAMIC_WINDOW_PARAMS.items()}
		dynamic_window_params.pop("dynamic_window")
		output_file = os.path.join(data_dir, 'actual_one_gene')
		plotfiles= [os.path.join(data_dir, 'valid'),os.path.join(data_dir, 'valid_ctrl')]
		pef= PrepareEMBLFile(plotfiles, 1, 4, 2, 15, os.path.join(data_dir, 'one_gene'),True,False,dynamic_window_params)
		plot_parser_objs = pef.plot_parser()
		eg = EMBLExpandGenes(os.path.join(data_dir, 'one_gene') , 15,True,False,dynamic_window_params)
		eg.construct_file(output_file, plot_parser_objs)
		self.assertTrue(os.path.exists(output_file))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_one_gene_dw'), output_file))
		os.remove(output_file)	

	#expand genes without dynamic window
	def test_embl_expand_genes(self):
		dynamic_window_params= {key: None for key, value in DYNAMIC_WINDOW_PARAMS.items()}
		dynamic_window_params.pop("dynamic_window")
		output_file = os.path.join(data_dir, 'actual_one_gene')
		plotfiles= [os.path.join(data_dir, 'valid'),os.path.join(data_dir, 'valid_ctrl')]
		pef= PrepareEMBLFile(plotfiles, 1, 4, 2, 15, os.path.join(data_dir, 'one_gene'),False,True,dynamic_window_params)
		plot_parser_objs = pef.plot_parser()
		e = EMBLExpandGenes(os.path.join(data_dir, 'one_gene') , 15,False,True,dynamic_window_params)
		e.construct_file(output_file, plot_parser_objs)
		self.assertTrue(os.path.exists(output_file))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_one_gene'), output_file))
		os.remove(output_file)		
		