# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for preparing reads
"""
import filecmp
import os
import unittest

from quatradis import mapper


class MapperTest(unittest.TestCase):
    """
    These unit tests are for the mapper functions
    """

    def test_calc_read_len(self):
        read_len = mapper.calc_read_length("data/mapper/test.fastq")
        self.assertEqual(44, read_len)

    def test_smalt_index_cmd(self):
        index_cmd, exitcode = mapper.index_reference("data/mapper/smallref.fa", "test.ref", 44, mapper="smalt", dry_run=True)
        self.assertEqual("smalt index -k 13 -s 4 test.ref data/mapper/smallref.fa > /dev/null 2>&1", index_cmd)
        self.assertEqual(0, exitcode)

    def test_smalt_align_cmd(self):
        align_cmd, exitcode = mapper.map_reads("data/mapper/test.fastq", "data/mapper/smallref.fa", "test.ref", "mapped.out", 44, mapper="smalt", dry_run=True)
        self.assertEqual("smalt map -x -n 1 test.ref data/mapper/test.fastq 1> mapped.out 2> align.stderr", align_cmd)
        self.assertEqual(0, exitcode)

    def test_bwa_index_cmd(self):
        index_cmd, exitcode = mapper.index_reference("data/mapper/smallref.fa", "test.ref", 44, mapper="bwa", dry_run=True)
        self.assertEqual("bwa index data/mapper/smallref.fa > /dev/null 2>&1", index_cmd)
        self.assertEqual(0, exitcode)

    def test_bwa_align_cmd(self):
        align_cmd, exitcode = mapper.map_reads("data/mapper/test.fastq", "data/mapper/smallref.fa", "test.ref", "mapped.out", 44, mapper="bwa", dry_run=True)
        self.assertEqual("bwa mem -k 13 -t 1 data/mapper/smallref.fa data/mapper/test.fastq 1> mapped.out 2> align.stderr", align_cmd)
        self.assertEqual(0, exitcode)

    # This test needs smalt installed
    def test_smalt(self):
        mapper.index_and_align("data/mapper/test.fastq", "data/mapper/smallref.fa", "test.ref", "mapped.out", mapper="smalt")
        os.system("grep -v ^@ mapped.out > mapped.nohead.out && grep -v ^@ data/mapper/expected.smalt.sam > expected.smalt.nohead.sam")
        self.assertTrue(filecmp.cmp("mapped.nohead.out", "expected.smalt.nohead.sam"))
        os.remove("test.ref.sma")
        os.remove("test.ref.smi")
        os.system("rm mapped.*")
        os.remove("expected.smalt.nohead.sam")
        os.remove("align.stderr")

    # This test need bwa installed
    def test_bwa(self):
        mapper.index_and_align("data/mapper/test.fastq", "data/mapper/smallref.fa", "test.ref", "mapped.out", mapper="bwa")
        os.system("grep -v ^@ mapped.out > mapped.nohead.out && grep -v ^@ data/mapper/expected.bwa.sam > expected.bwa.nohead.sam")
        self.assertTrue(filecmp.cmp("mapped.nohead.out", "expected.bwa.nohead.sam"))
        os.system("rm mapped.*")
        os.remove("align.stderr")
        os.remove("expected.bwa.nohead.sam")
        os.system("rm data/mapper/smallref.fa.*")

    def test_sam2bam(self):
        mapper.sam2bam("data/mapper/expected.bwa.sam", "nice.bam")
        self.assertTrue(os.path.exists("nice.bam"))
        self.assertTrue(os.path.exists("nice.bam.bai"))
        self.assertTrue(os.path.exists("nice.bam.bamcheck"))
        os.system("rm *nice.bam*")

    def test_sam2bam_multithread(self):
        mapper.sam2bam("data/mapper/expected.bwa.sam", "nice.bam", threads=2)
        self.assertTrue(os.path.exists("nice.bam"))
        self.assertTrue(os.path.exists("nice.bam.bai"))
        self.assertTrue(os.path.exists("nice.bam.bamcheck"))
        os.system("rm *nice.bam*")
