# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for analyse insert sites
"""
import unittest
from quatradis import isp_analyse

class AnalyseInsertSitesTest(unittest.TestCase):

    def test_get_cds_coords(self):
        #cds_coords = isp_analyse.get_cds_locations(embl_file=)
        self.assertTrue(True)

    def test_read_in_plot_file(self):
        data = isp_analyse.read_in_plot_file("data/isp_analyse/quatradis_out.plot.CP009273.1.insert_site_plot.gz")
        self.assertTrue(True)

    def test_get_insert_sites_from_plots(self):
        data = isp_analyse.get_insert_sites_from_plots(["data/isp_analyse/quatradis_out.plot.CP009273.1.insert_site_plot.gz"], False)
        self.assertTrue(True)
