# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for combining plot files
"""
import filecmp
import os.path
import unittest

from quatradis.tisp import combine
from tests.py.tisp import DATA_DIR
 
data_dir = os.path.join(DATA_DIR, "combine")

class CombinePlotsTest(unittest.TestCase):

    def tearDown(self):
        if os.path.exists("combined"):
            os.system("rm -r combined")


    def test_combine(self):
        plotfile = os.path.join(data_dir, "comb_sample.txt")
        combine.combine(plotfile)

        # Test expected files exist
        self.assertTrue(os.path.exists("combined/first.insert_site_plot.gz"))
        self.assertTrue(os.path.exists("combined/second.insert_site_plot.gz"))
        self.assertTrue(os.path.exists("comb_sample.stats"))

        # Now check they have the expected contents
        os.system("gunzip -c combined/first.insert_site_plot.gz > combined/first.test.plot")
        os.system("gunzip -c combined/second.insert_site_plot.gz > combined/second.test.plot")
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'first.expected.plot'), 'combined/first.test.plot'))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'second.expected.plot'), 'combined/second.test.plot'))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'comb_expected.stats'), 'comb_sample.stats'))

        os.remove("comb_sample.stats")

    def test_combine_gzipped(self):
        plotfile = os.path.join(data_dir, "zip_comb_list.txt")
        combine.combine(plotfile)

        self.assertTrue(os.path.exists("combined/zip_combined.insert_site_plot.gz"))
        self.assertTrue(os.path.exists("combined/tabix_sorted.insert_site_plot.gz"))
        self.assertTrue(os.path.exists("combined/tabix_sorted.insert_site_plot.gz.tbi"))

        os.system("gunzip -c combined/zip_combined.insert_site_plot.gz > combined/zip_combined.test.plot")
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, "zip_comb_exp.plot"), "combined/zip_combined.test.plot"))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, "zip_comb_exp.stats"), "zip_comb_list.stats"))

        os.remove("zip_comb_list.stats")

    def test_custom_directory(self):
        plotfile = os.path.join(data_dir, "comb_sample.txt")
        combine.combine(plotfile, "comb_test")

        self.assertTrue(os.path.exists("comb_test"))
        self.assertTrue(os.path.exists("comb_test/first.insert_site_plot.gz"))
        self.assertTrue(os.path.exists("comb_test/second.insert_site_plot.gz"))

        os.system("rm -r comb_test")
        os.remove("comb_sample.stats")
