#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir
DATA_DIR=tests/data/util/tags

# Sanity test help messages
echo -n "Checking 'tradis utils tags' help messages ... "
./tradis utils tags --help > /dev/null
./tradis utils tags add --help > /dev/null
./tradis utils tags check --help > /dev/null
./tradis utils tags filter --help > /dev/null
./tradis utils tags remove --help > /dev/null
echo "ok"

echo -n "Checking 'tradis utils tags add' ... "
./tradis utils tags add -o test_temp/test.tr.bam tests/data/tisp/create/small_multi_sequence.bam > /dev/null 2>&1 && rm -r test_temp
./tradis utils tags add tests/data/tisp/create/small_multi_sequence.bam > /dev/null 2>&1 && rm tests/data/tisp/create/small_multi_sequence.tr.bam
echo "ok"

echo -n "Checking 'tradis utils tags check' ... "
./tradis utils tags check $DATA_DIR/sample_sm_tr.bam 2>&1 | grep True > /dev/null
echo "ok"

echo -n "Checking 'tradis utils tags filter' ... "
./tradis utils tags filter --tag CAACGTTTT $DATA_DIR/sample.caa.fastq.gz test_temp/output.fastq.gz > /dev/null 2>&1 && rm -r test_temp
echo "ok"

echo -n "Checking 'tradis utils tags remove' ... "
./tradis utils tags remove --tag CAACGTTTT $DATA_DIR/sample.caa.fastq.gz test_temp/output.fastq.gz > /dev/null 2>&1 && rm -r test_temp
echo "ok"