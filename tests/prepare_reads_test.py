
# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for preparing reads
"""
import filecmp
import os
import unittest

from biotradis2 import prepare_reads


class PrepareReadsTest(unittest.TestCase):
    """
    These unit tests are for the prepare_reads functions
    """

    # ##### Tests for finding the tag

    def test_find_tag_exact(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTT")
        self.assertTrue(found)

    def test_find_tag_1mm(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTATT", max_mismatches=1)
        self.assertTrue(found)

    def test_find_notag_exact(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTATT", max_mismatches=0)
        self.assertFalse(found)

    def test_find_notag_1mm(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "TNAGAGACAG", max_mismatches=1)
        self.assertFalse(found)

    def test_find_lots_of_mismatches(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "TNAGAGACAG", max_mismatches=50)
        self.assertTrue(found)

    def test_find_long_tag(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTT")
        self.assertFalse(found)

    def test_find_long_tag2(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCACAACGTTTT")
        self.assertFalse(found)

    def test_find_notag(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "")
        self.assertFalse(found)

    def test_find_too_short_tag(self):
        found = prepare_reads.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "AT")
        self.assertFalse(found)

    # ##### Tests for file processing

    def test_file_exact_match(self):
        prepare_reads.prepare_reads("data/prepare_reads/sample.caa.fastq", "output.fastq", "CAACGTTTT")
        self.assertTrue(filecmp.cmp('output.fastq', 'data/prepare_reads/expected.rm.caa.fastq'))
        os.remove("output.fastq")

    def test_file_1mm(self):
        prepare_reads.prepare_reads("data/prepare_reads/sample.caa.fastq", "output.fastq", "CAACGTTTT", 1)
        self.assertTrue(filecmp.cmp('output.fastq', 'data/prepare_reads/expected.rm.1mm.caa.fastq'))
        os.remove("output.fastq")

    def test_file_tna(self):
        prepare_reads.prepare_reads("data/prepare_reads/sample.tna.fastq", "output.fastq", "TNAGAGACAG", 1)
        self.assertTrue(filecmp.cmp('output.fastq', 'data/prepare_reads/expected.rm.tna.fastq'))
        os.remove("output.fastq")

    def test_file_tna(self):
        try:
            prepare_reads.prepare_reads("data/prepare_reads/doesnotexist.fastq", "output.fastq", "TNAGAGACAG", 1)
            self.assertTrue(False)
        except FileNotFoundError:
            self.assertTrue(True)

    def test_file_exact_match_gzin(self):
        prepare_reads.prepare_reads("data/prepare_reads/sample.caa.fastq.gz", "output.fastq", "CAACGTTTT")
        self.assertTrue(filecmp.cmp('output.fastq', 'data/prepare_reads/expected.rm.caa.fastq'))
        os.remove("output.fastq")

    def test_file_exact_match_gzout(self):
        prepare_reads.prepare_reads("data/prepare_reads/sample.caa.fastq.gz", "output.fastq.gz", "CAACGTTTT")
        os.system("gunzip output.fastq.gz")
        os.system("gunzip -c data/prepare_reads/expected.rm.caa.fastq.gz > expected.fastq")
        self.assertTrue(filecmp.cmp('output.fastq', 'expected.fastq'))
        os.remove("output.fastq")
        os.remove("expected.fastq")


if __name__ == '__main__':
    unittest.main()