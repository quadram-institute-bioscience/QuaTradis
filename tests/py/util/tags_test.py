# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for general tradis tag handling, e.g. identification and removal in reads, and detection in alignments
"""
import filecmp
import os
import shutil
import unittest

import pysam
from pysam.libcalignedsegment import CIGAR_OPS

from quatradis.util import tags
from tests.py.util import DATA_DIR, CWD

data_dir = os.path.join(DATA_DIR, 'tags')
tisp_data_dir = os.path.join(CWD, '..', '..', 'data', 'tisp', 'create')


class TagsTest(unittest.TestCase):

    # ##### Tests for finding the tag

    def test_find_tag(self):
        tests = [
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTT", 0, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTATT", 1, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTATT", 0, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "TNAGAGACAG", 1, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "TNAGAGACAG", 50, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTTCAACGTTTT", 0, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCACAACGTTTT", 0, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "", 0, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", "AT", 0, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*CGT", 0, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*..CGT", 0, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGT", 0, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGTA", 0, False),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGT...CTG", 0, True),
            ("CAACGTTTTCTGCGTGTTGCCGATATTTTGGAAAGCA", ".*...CGT...CTGA", 0, False)

        ]
        for i, t in enumerate(tests):
            with self.subTest(i=i):
                found = tags.find_tag(t[0], t[1], max_mismatches=t[2])
                self.assertEqual(found, t[3])

    def test_remove_tags(self):
        tests = [
            ("sample.caa.fastq", "CAACGTTTT", 0, "expected.rm.caa.fastq", True),
            ("sample.caa.fastq", "CAACGTTTT", 1, "expected.rm.1mm.caa.fastq", True),
            ("sample.tna.fastq", "TNAGAGACAG", 1, "expected.rm.tna.fastq", True),
            ("doesnotexist.fastq", "TNAGAGACAG", 1, "expected.does.not.exist.fastq", False),
            ("sample.caa.fastq.gz", "CAACGTTTT", 0, "expected.rm.caa.fastq", True)
        ]
        for i, t in enumerate(tests):
            with self.subTest(i=i):
                try:
                    tags.remove_tags(os.path.join(data_dir, t[0]), "output.fastq", tag=t[1], max_mismatches=t[2])
                    self.assertTrue(t[4])
                    self.assertTrue(filecmp.cmp('output.fastq', os.path.join(data_dir, t[3])))
                    os.remove("output.fastq")
                except FileNotFoundError:
                    self.assertFalse(t[4])

    def test_file_exact_match_gzout(self):
        tags.remove_tags(os.path.join(data_dir, "sample.caa.fastq.gz"), "output.fastq.gz", "CAACGTTTT")
        os.system("gunzip output.fastq.gz")
        os.system("gunzip -c " + os.path.join(data_dir, "expected.rm.caa.fastq.gz") + " > expected.fastq")
        self.assertTrue(filecmp.cmp('output.fastq', 'expected.fastq'))
        os.remove("output.fastq")
        os.remove("expected.fastq")

    def testTagsInAlignment(self):
        tests = [
            ("sample_sm_tr.bam", True),
            ("sample_sm_no_tr.bam", False)
            #("sample_sm_tr.cram", True) # This won't work until we regenerate a suitable cram file
        ]
        for i, t in enumerate(tests):
            with self.subTest(i=i):
                bamfile = os.path.join(data_dir, t[0])
                res = tags.tags_in_alignment(bamfile)
                self.assertEqual(res, t[1])

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
        bamfile = os.path.join(tisp_data_dir, "small_multi_sequence.bam")
        outfile = "test_temp/test.tr.bam"
        tags.add_tags(bamfile, outfile)
        self.assertTrue(True)
        shutil.rmtree("test_temp")

    def testAddTagsToBamNoOutfile(self):
        bamfile=os.path.join(tisp_data_dir, "small_multi_sequence.bam")
        tags.add_tags(bamfile)
        self.assertTrue(True)
        os.remove(os.path.join(tisp_data_dir, "small_multi_sequence.tr.bam"))


if __name__ == '__main__':
    unittest.main()