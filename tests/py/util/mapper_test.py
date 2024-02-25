# !/usr/bin/env python3
# coding: utf-8

"""
Unit tests for preparing reads
"""
import filecmp
import os
import unittest

from quatradis.util import mapper

data_dir = os.path.join("data", "util", "mapper")


class MapperTest(unittest.TestCase):
    """
    These unit tests are for the mapper functions
    """

    def test_calc_read_len(self):
        read_len = mapper.calc_read_length(os.path.join(data_dir, "test.fastq"))
        self.assertEqual(44, read_len)

    def test_smalt_index_cmd(self):
        index_cmd, exitcode = mapper.index_reference(
            os.path.join(data_dir, "smallref.fa"),
            "test.ref",
            44,
            mapper="smalt",
            dry_run=True,
        )
        self.assertEqual(
            "smalt index -k 13 -s 4 test.ref "
            + os.path.join(data_dir, "smallref.fa")
            + " > /dev/null 2>&1",
            index_cmd,
        )
        self.assertEqual(0, exitcode)

    def test_smalt_align_cmd(self):
        align_cmd, exitcode = mapper.map_reads(
            os.path.join(data_dir, "test.fastq"),
            os.path.join(data_dir, "smallref.fa"),
            "test.ref",
            "mapped.out",
            44,
            mapper="smalt",
            dry_run=True,
        )
        self.assertEqual(
            "smalt map -x -n 1 test.ref "
            + os.path.join(data_dir, "test.fastq")
            + " 1> mapped.out 2> mapped.out.stderr",
            align_cmd,
        )
        self.assertEqual(0, exitcode)

    def test_bwa_index_cmd(self):
        index_cmd, exitcode = mapper.index_reference(
            os.path.join(data_dir, "smallref.fa"),
            "test.ref",
            44,
            mapper="bwa",
            dry_run=True,
        )
        self.assertEqual(
            "bwa index " + os.path.join(data_dir, "smallref.fa") + " > /dev/null 2>&1",
            index_cmd,
        )
        self.assertEqual(0, exitcode)

    def test_bwa_align_cmd(self):
        align_cmd, exitcode = mapper.map_reads(
            os.path.join(data_dir, "test.fastq"),
            os.path.join(data_dir, "smallref.fa"),
            "test.ref",
            "mapped.out",
            44,
            mapper="bwa",
            dry_run=True,
        )
        self.assertEqual(
            "bwa mem -k 13 -t 1 "
            + os.path.join(data_dir, "smallref.fa")
            + " "
            + os.path.join(data_dir, "test.fastq")
            + " 1> mapped.out 2> mapped.out.stderr",
            align_cmd,
        )
        self.assertEqual(0, exitcode)

    # This test needs smalt installed
    def test_smalt(self):
        mapper.index_and_align(
            os.path.join(data_dir, "test.fastq"),
            os.path.join(data_dir, "smallref.fa"),
            "test.ref",
            "mapped.out",
            mapper="smalt",
        )
        os.system(
            "grep -v ^@ mapped.out > mapped.nohead.out && grep -v ^@ "
            + os.path.join(data_dir, "expected.smalt.sam")
            + " > expected.smalt.nohead.sam"
        )
        self.assertTrue(filecmp.cmp("mapped.nohead.out", "expected.smalt.nohead.sam"))
        os.remove("test.ref.sma")
        os.remove("test.ref.smi")
        os.system("rm mapped.*")
        os.remove("expected.smalt.nohead.sam")

    # This test need bwa installed
    def test_bwa(self):
        mapper.index_and_align(
            os.path.join(data_dir, "test.fastq"),
            os.path.join(data_dir, "smallref.fa"),
            "test.ref",
            "mapped.out",
            mapper="bwa",
        )
        os.system(
            "grep -v ^@ mapped.out > mapped.nohead.out && grep -v ^@ "
            + os.path.join(data_dir, "expected.bwa.sam")
            + " > expected.bwa.nohead.sam"
        )
        self.assertTrue(filecmp.cmp("mapped.nohead.out", "expected.bwa.nohead.sam"))
        os.system("rm mapped.*")
        os.remove("expected.bwa.nohead.sam")
        os.system("rm " + os.path.join(data_dir, "smallref.fa") + ".*")

    def test_sam2bam(self):
        mapper.sam2bam(os.path.join(data_dir, "expected.bwa.sam"), "nice.bam")
        self.assertTrue(os.path.exists("nice.bam"))
        self.assertTrue(os.path.exists("nice.bam.bai"))
        self.assertTrue(os.path.exists("nice.bam.bamcheck"))
        os.system("rm *nice.bam*")

    def test_sam2bam_multithread(self):
        mapper.sam2bam(
            os.path.join(data_dir, "expected.bwa.sam"), "nice.bam", threads=2
        )
        self.assertTrue(os.path.exists("nice.bam"))
        self.assertTrue(os.path.exists("nice.bam.bai"))
        self.assertTrue(os.path.exists("nice.bam.bamcheck"))
        os.system("rm *nice.bam*")
