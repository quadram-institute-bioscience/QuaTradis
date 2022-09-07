import unittest

from quatradis.embl.window import WindowGenerator


class TestWindowGenerator(unittest.TestCase):

    def test_window_generator_nonoverlap(self):
        w = WindowGenerator(20, 5, 5)
        self.assertEqual(4, len(w.create_windows()))

    def test_window_generator_overlap(self):
        w = WindowGenerator(20, 5, 1)
        self.assertEqual(16, len(w.create_windows()))
