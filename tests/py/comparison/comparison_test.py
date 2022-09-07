import filecmp
import os
import shutil
import unittest

from quatradis.comparison.comparison import prep_essentiality_pipeline_output_for_comparison, run_comparisons
from quatradis.comparison.plot_logfc import PlotLogOptions

data_dir = os.path.join('data', 'essentiality', 'output')
data_comp_dir = os.path.join('data', 'comparison')



class TestComparison(unittest.TestCase):

    def setUp(self):
        self.embl_file = os.path.join(data_dir, 'prepared.embl')



    def test_prep(self):
        controls_all, controls_only_ess, conditions_all, conditions_only_ess = prep_essentiality_pipeline_output_for_comparison(
            [os.path.join(data_dir, 'small_control')], [os.path.join(data_dir, 'small_case')], 'combined')

        self.assertEquals(['data/essentiality/output/small_case/combined.plot.gz.stats.all.tsv'], conditions_all)
        self.assertEquals(['data/essentiality/output/small_control/combined.plot.gz.stats.all.tsv'], controls_all)
        self.assertEquals(['data/essentiality/output/small_case/combined.plot.gz.stats.essen.csv'], conditions_only_ess)
        self.assertEquals(['data/essentiality/output/small_control/combined.plot.gz.stats.essen.csv'], controls_only_ess)

    def test_small_real(self):
        controls_all, controls_only_ess, conditions_all, conditions_only_ess = prep_essentiality_pipeline_output_for_comparison(
            [os.path.join(data_dir, 'small_control')], [os.path.join(data_dir, 'small_case')], 'combined')

        options = PlotLogOptions(genome_length=50000, minimum_logfc=1, maximum_pvalue=1, maximum_qvalue=1, minimum_logcpm=1, window_size=100,
                       span_gaps=0, report_decreased_insertions=True, minimum_block=1)

        run_comparisons(controls_all, conditions_all, controls_only_ess, conditions_only_ess, self.embl_file, comparison_options=options, prefix='testoutput', verbose=False)

        self.assertTrue(os.path.exists('testoutput'))
        shutil.rmtree("testoutput")

    def test_ignore_decreased_insertions(self):
        controls_all, controls_only_ess, conditions_all, conditions_only_ess = prep_essentiality_pipeline_output_for_comparison(
            [os.path.join(data_dir, 'small_control_high_insertions')], [os.path.join(data_dir, 'small_case')], 'combined')

        options = PlotLogOptions(genome_length=50000, minimum_logfc=1, maximum_pvalue=1, maximum_qvalue=1,
                                 minimum_logcpm=1, window_size=100,
                                 span_gaps=0, report_decreased_insertions=False, minimum_block=1)

        run_comparisons(controls_all, conditions_all, controls_only_ess, conditions_only_ess, self.embl_file,
                        comparison_options=options, prefix='testoutputx', verbose=False)

        self.assertTrue(os.path.exists('testoutputx'))
        shutil.rmtree("testoutputx")


