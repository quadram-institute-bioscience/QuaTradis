# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for analyse insert sites
"""
import os
import unittest

from Bio.SeqFeature import SeqFeature

from quatradis import isp_analyse

class AnalyseInsertSitesTest(unittest.TestCase):

    def test_get_cds_coords(self):
        cds_coords = isp_analyse.get_cds_locations(embl_file="data/isp_analyse/reference_BW25113_short.embl")
        self.assertEqual(922, len(cds_coords))

    def test_read_in_plot_file(self):
        data = isp_analyse.read_in_plot_file("data/isp_analyse/controlLBrep1.insert_site_plot_short.gz")
        self.assertEqual(1000002, len(data))

    def test_get_insert_sites_from_plots(self):
        data = isp_analyse.get_insert_sites_from_plots(["data/isp_analyse/controlLBrep1.insert_site_plot_short.gz"], False)
        self.assertEqual(1, len(data))
        self.assertEqual(2, data.ndim)
        self.assertEqual(1000002, len(data[0]))

    def test_analyse_insert_sites(self):
        isp_analyse.analyse_insert_sites(embl_file="data/isp_analyse/reference_BW25113_short.embl",
                                                plot_files=["data/isp_analyse/controlLBrep1.insert_site_plot_short.gz",
                                                            "data/isp_analyse/025mgLTricRep1.insert_site_plot_short.gz"])
        self.assertTrue(os.path.exists("controlLBrep1.insert_site_plot_short.gz.tradis_gene_insert_sites.csv"))
        self.assertTrue(os.path.exists("025mgLTricRep1.insert_site_plot_short.gz.tradis_gene_insert_sites.csv"))
        os.system("rm *.tradis_gene_insert_sites.csv")

    def test_analyse_insert_sites_joined(self):
        isp_analyse.analyse_insert_sites(embl_file="data/isp_analyse/reference_BW25113_short.embl",
                                                plot_files=["data/isp_analyse/controlLBrep1.insert_site_plot_short.gz",
                                                            "data/isp_analyse/025mgLTricRep1.insert_site_plot_short.gz"],
                                         joined_output=True)
        self.assertTrue(os.path.exists("joined_output.tradis_gene_insert_sites.csv"))
        os.system("rm *.tradis_gene_insert_sites.csv")

    def test_gene_name_good(self):
        f = SeqFeature()
        f.qualifiers['gene'] = ["imagene"]
        a = isp_analyse.get_gene_name(f)
        self.assertEqual('imagene', a)

    def test_gene_name_withsymbols(self):
        f = SeqFeature()
        f.qualifiers['gene'] = ["image$ [{&!  ne"]
        a = isp_analyse.get_gene_name(f)
        self.assertEqual('imagene', a)



