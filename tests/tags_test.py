# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for general tradis tag handling, e.g. identification and removal in reads, and detection in alignments
"""
import filecmp
import os
import unittest

from quatradis import tags


class TagsTest(unittest.TestCase):

    # ##### Tests for finding the tag

    def test_find_tag_exact(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTT")
        self.assertTrue(found)

    def test_find_tag_1mm(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTATT", max_mismatches=1)
        self.assertTrue(found)

    def test_find_notag_exact(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTATT", max_mismatches=0)
        self.assertFalse(found)

    def test_find_notag_1mm(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "TNAGAGACAG", max_mismatches=1)
        self.assertFalse(found)

    def test_find_lots_of_mismatches(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "TNAGAGACAG", max_mismatches=50)
        self.assertTrue(found)

    def test_find_long_tag(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTT")
        self.assertFalse(found)

    def test_find_long_tag2(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCACAACGTTTT")
        self.assertFalse(found)

    def test_find_notag(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "")
        self.assertFalse(found)

    def test_find_too_short_tag(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "AT")
        self.assertFalse(found)

    def test_find_tag_regex_1(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*CGT")
        self.assertTrue(found)

    def test_find_tag_regex_2(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*..CGT")
        self.assertTrue(found)

    def test_find_tag_regex_3(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGT")
        self.assertTrue(found)

    def test_find_tag_regex_4(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGTA")
        self.assertFalse(found)

    def test_find_tag_regex_5(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGT...CTG")
        self.assertTrue(found)

    def test_find_tag_regex_6(self):
        found = tags.find_tag("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGT...CTGA")
        self.assertFalse(found)

    # ##### Tests for file processing

    def test_file_exact_match(self):
        tags.remove_tags("data/tags/sample.caa.fastq", "output.fastq", "CAACGTTTT")
        self.assertTrue(filecmp.cmp('output.fastq', 'data/tags/expected.rm.caa.fastq'))
        os.remove("output.fastq")

    def test_file_1mm(self):
        tags.remove_tags("data/tags/sample.caa.fastq", "output.fastq", tag="CAACGTTTT", max_mismatches=1)
        self.assertTrue(filecmp.cmp('output.fastq', 'data/tags/expected.rm.1mm.caa.fastq'))
        os.remove("output.fastq")

    def test_file_tna(self):
        tags.remove_tags("data/tags/sample.tna.fastq", "output.fastq", tag="TNAGAGACAG", max_mismatches=1)
        self.assertTrue(filecmp.cmp('output.fastq', 'data/tags/expected.rm.tna.fastq'))
        os.remove("output.fastq")

    def test_file_tna(self):
        try:
            tags.remove_tags("data/tags/doesnotexist.fastq", "output.fastq", tag="TNAGAGACAG", max_mismatches=1)
            self.assertTrue(False)
        except FileNotFoundError:
            self.assertTrue(True)

    def test_file_exact_match_gzin(self):
        tags.remove_tags("data/tags/sample.caa.fastq.gz", "output.fastq", "CAACGTTTT")
        self.assertTrue(filecmp.cmp('output.fastq', 'data/tags/expected.rm.caa.fastq'))
        os.remove("output.fastq")

    def test_file_exact_match_gzout(self):
        tags.remove_tags("data/tags/sample.caa.fastq.gz", "output.fastq.gz", "CAACGTTTT")
        os.system("gunzip output.fastq.gz")
        os.system("gunzip -c data/tags/expected.rm.caa.fastq.gz > expected.fastq")
        self.assertTrue(filecmp.cmp('output.fastq', 'expected.fastq'))
        os.remove("output.fastq")
        os.remove("expected.fastq")


    # Test tag detection in SAM/BAM/CRAM

    def testGoodBam(self):
        bamfile = "data/tags/sample_sm_tr.bam"
        self.assertTrue(tags.tags_in_alignment(bamfile))

    # This won't work until we regenerate a suitable cram file
    # def testGoodCram(self):
    #     cramfile = "data/tags/sample_sm_tr.cram"
    #     self.assertTrue(tags.tags_in_alignment(cramfile))

    def testNoTradis(self):
        bamfile = "data/tags/sample_sm_no_tr.bam"
        self.assertFalse(tags.tags_in_alignment(bamfile))

if __name__ == '__main__':
    unittest.main()