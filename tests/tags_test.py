# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for general tradis tag handling, e.g. identification and removal in reads, and detection in alignments
"""
import filecmp
import os
import unittest

import pysam
from pysam.libcalignedsegment import CIGAR_OPS

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


    # Add tags

    def testAddTag(self):
        a = pysam.AlignedSegment()
        a.query_name = "read_28833_29006_6945"
        a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        a.flag = 99
        a.reference_id = 0
        a.reference_start = 32
        a.mapping_quality = 20
        a.cigar = ((0, 10), (2, 1), (0, 25))
        a.next_reference_id = 0
        a.next_reference_start = 199
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
        a.tags = (("NM", 1),
                  ("RG", "L1"),
                  ("tr", "GGGG"),
                  ("tq", "8888"))

        tagged = tags.add_tags_to_alignment(a)
        self.assertEqual((CIGAR_OPS.CMATCH.value, 39), tagged.cigartuples[0])
        self.assertEqual(1, len(tagged.cigartuples))
        self.assertEqual('8888<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<', tagged.qual)
        self.assertEqual('GGGGAGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG', tagged.query_sequence)

    def testAddTagToReversed(self):
        a = pysam.AlignedSegment()
        a.query_name = "read_28833_29006_6945"
        a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
        a.flag = 99
        a.reference_id = 0
        a.reference_start = 32
        a.is_reverse = True
        a.mapping_quality = 20
        a.cigar = ((0, 10), (2, 1), (0, 25))
        a.next_reference_id = 0
        a.next_reference_start = 199
        a.template_length = 167
        a.tags = (("NM", 1),
                  ("RG", "L1"),
                  ("tr", "GGAA"),
                  ("tq", "8877"))

        tagged = tags.add_tags_to_alignment(a)
        self.assertEqual((CIGAR_OPS.CMATCH.value, 39), tagged.cigartuples[0])
        self.assertEqual(1, len(tagged.cigartuples))
        self.assertEqual('AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCGTTCC', tagged.query_sequence)
        self.assertEqual('<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<7788', tagged.qual)

    def testAddTagsToBam(self):
        bamfile="data/isp_create/small_multi_sequence.bam"
        outfile="test.tr.bam"
        tags.add_tags(bamfile, outfile)
        self.assertTrue(True)
        os.remove(outfile)

    def testAddTagsToBamNoOutfile(self):
        bamfile="data/isp_create/small_multi_sequence.bam"
        tags.add_tags(bamfile)
        self.assertTrue(True)
        os.remove("data/isp_create/small_multi_sequence.tr.bam")

if __name__ == '__main__':
    unittest.main()