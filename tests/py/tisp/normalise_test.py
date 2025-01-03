import filecmp
import os
import unittest

import shutil

from quatradis.tisp.normalise import NormalisePlots
from tests.py.tisp import DATA_DIR

data_dir = os.path.join(DATA_DIR, "normalise")
plotparser_dir = os.path.join(DATA_DIR, "parser")


class NormaliseTest(unittest.TestCase):

    def test_big_differences(self):
        n = NormalisePlots([os.path.join(data_dir, "sample1"), os.path.join(data_dir, "sample2")],
                           0.1, output_dir="normalised")
        output_files, max_reads = n.create_normalised_files()
        self.assertEqual(2, len(output_files))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, "sample2"), output_files[1]))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, "expected_sample1"), output_files[0]))
        self.assertTrue(n.decreased_insertion_reporting())
        shutil.rmtree("normalised")

    def test_decreased_insertion_reporting(self):
        tests = [
            ([os.path.join(data_dir, "lowinsertions"), os.path.join(data_dir, "highinsertions")], False),
            ([os.path.join(data_dir, "fewinsertions"), os.path.join(data_dir, "manyinsertions")], False),
            ([os.path.join(plotparser_dir, "Control2.out.CP009273.insert_site_plot.gz"),
              os.path.join(plotparser_dir, "Chloramrep2-MICpool.out.CP009273.insert_site_plot.gz")], True)
        ]

        for i, t in enumerate(tests):
            with self.subTest(i=i):
                n = NormalisePlots(t[0], 0.1, output_dir="normalised")
                if t[1]:
                    output_files, max_reads = n.create_normalised_files()
                actual = n.decreased_insertion_reporting()
                self.assertEqual(actual, t[1])
                if actual:
                    shutil.rmtree("normalised")

