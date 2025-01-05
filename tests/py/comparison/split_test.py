import filecmp
import shutil
import unittest
import os
from quatradis.comparison.split import split_plot
from tests.py.comparison import DATA_DIR

data_dir = os.path.join(DATA_DIR, 'split')


class TestSplit(unittest.TestCase):

    def test_split_plot(self):
        output_dir = os.path.join(data_dir, "small_case", "output")

        split_plot(os.path.join(data_dir, 'small_case.insert_site_plot.gz'), output_dir, minimum_threshold=5)

        self.assertTrue(os.path.exists(os.path.join(output_dir, "combined.plot.gz")))
        self.assertTrue(filecmp.cmp(os.path.join(output_dir, "combined.plot.gz"),
                                    os.path.join(data_dir, "small_case", "combined.plot.gz")))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "forward.plot.gz")))
        self.assertTrue(filecmp.cmp(os.path.join(output_dir, "forward.plot.gz"),
                                    os.path.join(data_dir, "small_case", "forward.plot.gz")))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "reverse.plot.gz")))
        self.assertTrue(filecmp.cmp(os.path.join(output_dir, "reverse.plot.gz"),
                                    os.path.join(data_dir, "small_case", "reverse.plot.gz")))

        shutil.rmtree(output_dir)
