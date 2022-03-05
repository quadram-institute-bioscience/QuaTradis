#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir

# Sanity test help messages
echo -n "Checking 'tradis tag' help messages ... "
./scripts/tradis tags --help > /dev/null
./scripts/tradis tags add --help > /dev/null
./scripts/tradis tags check --help > /dev/null
./scripts/tradis tags filter --help > /dev/null
./scripts/tradis tags remove --help > /dev/null
echo "ok"

echo -n "Checking 'tradis tag add' ... "
./scripts/tradis tags add -o test_temp/test.tr.bam tests/data/isp_create/small_multi_sequence.bam 2>&1 > /dev/null && rm -r test_temp
./scripts/tradis tags add tests/data/isp_create/small_multi_sequence.bam 2>&1 > /dev/null && rm tests/data/isp_create/small_multi_sequence.tr.bam
echo "ok"

echo -n "Checking 'tradis tag check' ... "
./scripts/tradis tags check tests/data/tags/sample_sm_tr.bam 2>&1 | grep True > /dev/null
echo "ok"

echo -n "Checking 'tradis tag filter' ... "
./scripts/tradis tags filter --tag CAACGTTTT tests/data/tags/sample.caa.fastq.gz test_temp/output.fastq.gz 2>&1 > /dev/null && rm -r test_temp
echo "ok"

echo -n "Checking 'tradis tag remove' ... "
./scripts/tradis tags remove --tag CAACGTTTT tests/data/tags/sample.caa.fastq.gz test_temp/output.fastq.gz 2>&1 > /dev/null && rm -r test_temp
echo "ok"