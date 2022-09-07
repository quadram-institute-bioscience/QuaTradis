import filecmp
import os
import unittest

from quatradis.tisp.normalise import NormalisePlots


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass


data_dir = os.path.join('data', 'tisp', 'normalise')
plotparser_dir = os.path.join('data', 'tisp', 'parser')


class TestNormalisePlots(unittest.TestCase):

    def test_big_differences(self):
        p = NormalisePlots([os.path.join(data_dir, 'sample1'), os.path.join(data_dir, 'sample2')], 0.1)
        output_files, max_reads = p.create_normalised_files()
        self.assertEqual(2, len(output_files))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sample2'), output_files[1]))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_sample1'), output_files[0]))
        self.assertTrue(p.decreased_insertion_reporting())

    def test_ignore_decreased_insertions(self):
        p = NormalisePlots([os.path.join(data_dir, 'lowinsertions'), os.path.join(data_dir, 'highinsertions')], 0.1)
        self.assertFalse(p.decreased_insertion_reporting())

    def test_ignore_decreased_insertions_insert_sites(self):
        p = NormalisePlots([os.path.join(data_dir, 'fewinsertions'), os.path.join(data_dir, 'manyinsertions')], 0.1)
        self.assertFalse(p.decreased_insertion_reporting())

    def test_large_file_zipped(self):
        p = NormalisePlots([os.path.join(plotparser_dir, 'Control2.out.CP009273.insert_site_plot.gz'),
                            os.path.join(plotparser_dir, 'Chloramrep2-MICpool.out.CP009273.insert_site_plot.gz')], 0.1)
        output_files, max_reads = p.create_normalised_files()
        self.assertTrue(p.decreased_insertion_reporting())
