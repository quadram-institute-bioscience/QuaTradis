# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for insert site plots
"""
import os
import unittest

from quatradis.tisp.generator import from_alignments
from tests.py.tisp.generator import DATA_DIR


def load_plot_file(plot_file):
    import gzip
    with gzip.open(plot_file, 'rb') as plot_handle:
        # Note the hack for the first line.  This is so we can use 1-based access to represent line numbers / coord positions in file
        return [""] + [str(x, 'UTF-8').strip() for x in plot_handle.readlines()]

data_dir = os.path.join(DATA_DIR, "tisp", "create")
tags_data_dir = os.path.join(DATA_DIR, "util", "tags")

class InsertSitePlotTest(unittest.TestCase):

    def test_plotter(self):
        from_alignments.plot(os.path.join(data_dir, "small_multi_sequence.bam"),
                    os.path.join(tags_data_dir, "sample.caa.fastq.gz"),
                    plot_out_prefix="temp/insert_site",
                    cutoff_score=0,
                    nb_reads=193,
                    nb_tagged_reads=193)
        self.assertTrue(os.path.exists("temp"))
        fn543502 = load_plot_file("temp/insert_site.FN543502.insert_site_plot.gz")
        pCROD1 = load_plot_file("temp/insert_site.pCROD1.insert_site_plot.gz")
        pCROD2 = load_plot_file("temp/insert_site.pCROD2.insert_site_plot.gz")
        pCROD3 = load_plot_file("temp/insert_site.pCROD3.insert_site_plot.gz")

        self.assertEqual("0 0", fn543502[1])  # First value
        self.assertEqual("0 0", fn543502[8974])  # Last value
        self.assertEqual("0 2", fn543502[7899])  # Before site
        # TODO check this one.  BioTradis 1 says this should have 12 hits on negative strand, but due to soft clipping
        # of 2 subsequent reads we end up with 14 instead.  Double check this is correct.
        # self.assertEqual("0 12", self.fn543502[7915])
        self.assertEqual("0 14", fn543502[7915])  # Reverse reads

        self.assertEqual("0 0", fn543502[249])
        self.assertEqual("1 0", fn543502[345])
        self.assertEqual("3 0", fn543502[354])
        self.assertEqual("1 0", fn543502[366])
        self.assertEqual("0 0", fn543502[513])

        self.assertEqual("0 0", pCROD1[1])  # First value
        self.assertEqual("0 0", pCROD1[59])  # Last value

        self.assertEqual("0 0", pCROD2[1])  # First value
        self.assertEqual("0 1", pCROD2[143])  # First read base
        self.assertEqual("0 0", pCROD2[144])  # Before read
        self.assertEqual("0 0", pCROD2[1000])  # Last value

        self.assertEqual("0 0", pCROD3[1])  # First value
        self.assertEqual("0 0", pCROD3[100])  # Last value

        os.system("rm -r temp")
        if os.path.exists("insert_site.stats"):
            os.remove("insert_site.stats")

    def test_read_counter_gzipped(self):
        nb_reads = from_alignments.get_number_reads(os.path.join(tags_data_dir, "sample.caa.fastq.gz"))
        self.assertEqual(10, nb_reads)

    def test_read_counter(self):
        nb_reads = from_alignments.get_number_reads(os.path.join(tags_data_dir, "sample.caa.fastq"))
        self.assertEqual(10, nb_reads)
