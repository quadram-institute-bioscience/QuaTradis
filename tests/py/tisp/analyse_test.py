# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for analyse insert sites
"""
import os
import unittest

from Bio.SeqFeature import SeqFeature

from quatradis.tisp import analyse
from tests.py.tisp import DATA_DIR

data_dir = os.path.join(DATA_DIR, "analyse")


class AnalyseInsertSitesTest(unittest.TestCase):

    def test_get_cds_coords(self):
        cds_coords = analyse.get_cds_locations(embl_file=os.path.join(data_dir, "reference_BW25113_short.embl"))
        self.assertEqual(922, len(cds_coords))

    def test_read_in_plot_file(self):
        data = analyse.read_in_plot_file(os.path.join(data_dir, "controlLBrep1.insert_site_plot_short.gz"))
        self.assertEqual(1000002, len(data))

    def test_get_insert_sites_from_plots(self):
        data = analyse.get_insert_sites_from_plots([os.path.join(data_dir, "controlLBrep1.insert_site_plot_short.gz")], False)
        self.assertEqual(1, len(data))
        self.assertEqual(2, data.ndim)
        self.assertEqual(1000002, len(data[0]))

    def test_analyse_insert_sites(self):
        analyse.count_insert_sites(embl_file=os.path.join(data_dir, "reference_BW25113_short.embl"),
                                   plot_files=[os.path.join(data_dir, "controlLBrep1.insert_site_plot_short.gz"),
                                                            os.path.join(data_dir, "025mgLTricRep1.insert_site_plot_short.gz")])
        self.assertTrue(os.path.exists("controlLBrep1.tradis_gene_insert_sites.csv"))
        self.assertTrue(os.path.exists("025mgLTricRep1.tradis_gene_insert_sites.csv"))
        os.system("rm *.tradis_gene_insert_sites.csv")

    def test_analyse_insert_sites_joined(self):
        analyse.count_insert_sites(embl_file=os.path.join(data_dir, "reference_BW25113_short.embl"),
                                   plot_files=[os.path.join(data_dir, "controlLBrep1.insert_site_plot_short.gz"),
                                                            os.path.join(data_dir, "025mgLTricRep1.insert_site_plot_short.gz")],
                                   joined_output=True)
        self.assertTrue(os.path.exists("joined_output.tradis_gene_insert_sites.csv"))
        os.system("rm *.tradis_gene_insert_sites.csv")

    def test_gene_name_good(self):
        f = SeqFeature()
        f.qualifiers['gene'] = ["imagene"]
        a = analyse.get_gene_name(f)
        self.assertEqual('imagene', a)

    def test_gene_name_withsymbols(self):
        f = SeqFeature()
        f.qualifiers['gene'] = ["image$ [{&!  ne"]
        a = analyse.get_gene_name(f)
        self.assertEqual('imagene', a)



