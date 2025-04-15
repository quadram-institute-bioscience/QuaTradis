import unittest
import os
import filecmp
from quatradis.comparison.gene_stats import gene_statistics_new, gene_statistics_old
from quatradis.util.config_defaults import GENE_REPORT_PARAMS
from quatradis.comparison.gene_stats import write_gene_report_new
from tests.py.comparison import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'gene_stats')

class TestGeneStats(unittest.TestCase):

    def test_new_gene_algo(self):
        # prepare input data
        plot_dir = os.path.join(data_dir, 'plot_files')
        all_files = sorted(os.listdir(plot_dir))  # Sorting ensures order

        # Get absolute paths and separate by condition
        cndtn_files = [os.path.join(plot_dir, f) for f in all_files if "cndtn" in f]
        ctrl_files = [os.path.join(plot_dir, f) for f in all_files if "ctrl" in f]
        plotfiles_all = cndtn_files + ctrl_files

        # Get absolute paths for all relevant files
        forward_count_condition = [os.path.join(data_dir, 'forward_count_condition_tsv', f) 
                                for f in os.listdir(os.path.join(data_dir, 'forward_count_condition_tsv'))]
        
        reverse_count_condition = [os.path.join(data_dir, 'reverse_count_condition_tsv', f) 
                                for f in os.listdir(os.path.join(data_dir, 'reverse_count_condition_tsv'))]
        
        forward_count_control = [os.path.join(data_dir, 'forward_count_control_tsv', f) 
                                for f in os.listdir(os.path.join(data_dir, 'forward_count_control_tsv'))]
        
        reverse_count_control = [os.path.join(data_dir, 'reverse_count_control_tsv', f) 
                                for f in os.listdir(os.path.join(data_dir, 'reverse_count_control_tsv'))]

        combined_compare = os.path.join(data_dir, 'compare_csv', 'combined.compare.csv')
        forward_compare = os.path.join(data_dir, 'compare_csv', 'forward.compare.csv')
        reverse_compare = os.path.join(data_dir, 'compare_csv', 'reverse.compare.csv')
        embl = os.path.join(data_dir, 'prepared.embl')
        output_file = os.path.join(data_dir,"gene_report.tsv")  # Already absolute

        gene_categorization_params_values = {key: value[2] for key, value in GENE_REPORT_PARAMS.items()}
        print("Gene Report Categorization Parameters:", gene_categorization_params_values)

        
        genereport= gene_statistics_new(
            False,
            plotfiles_all=plotfiles_all,
            forward_count_condition=forward_count_condition,
            reverse_count_condition=reverse_count_condition,
            forward_count_control=forward_count_control,
            reverse_count_control=reverse_count_control,
            combined_compare_csv=combined_compare,
            forward_compare_csv=forward_compare,
            reverse_compare_csv=reverse_compare,
            embl_file=embl,
            output_dir=output_file,
            gene_categorization_params_values=gene_categorization_params_values,
        )
        print("Check 0")
        #Assert Knockout
        self.assertEqual(3, (genereport['Category1'] == 'knockout').sum())
        print("Check 1")
        #Assert Protection
        self.assertEqual(7, (genereport['Category1'] == 'protection').sum())
        print("Check 2")
        #Assert Fractional Knockout
        self.assertEqual(7, (genereport['Category1'] == 'fractional knockout').sum())
        print("Check 3")
        #Assert Upregulation
        self.assertEqual(6, (genereport['Category2'] == 'upregulated').sum())
        print("Check 4")
        #Assert Downregulation
        self.assertEqual(2, (genereport['Category3'] == 'downregulated').sum())
        print("Check 5")

