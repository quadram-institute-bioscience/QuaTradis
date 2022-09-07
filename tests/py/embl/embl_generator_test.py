import unittest
import os

from quatradis.embl.generator import EMBLGenerator
from quatradis.embl.window import WindowGenerator

data_dir = os.path.join('data', 'embl', 'generator')


class TestEMBLGenerator(unittest.TestCase):

    def test_embl_generator_nonoverlap(self):
        w = WindowGenerator(20, 5, 5)
        e = EMBLGenerator(w.create_windows(), 20)

        e.construct_file(os.path.join(data_dir, 'nonoverlap'))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'nonoverlap')))
