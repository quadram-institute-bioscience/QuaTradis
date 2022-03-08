#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir

# Sanity test help messages
echo -n "Checking 'tradis pipeline' help messages ... "
./scripts/tradis pipeline --help > /dev/null
./scripts/tradis pipeline single --help > /dev/null
./scripts/tradis pipeline multiple --help > /dev/null
echo "ok"

echo -n "Checking 'tradis pipeline single' ... "
./scripts/tradis pipeline single --output_dir temp_test --profile tests/data/mapper/test.fastq tests/data/mapper/smallref.fa 2>&1 > /dev/null && rm tests/data/mapper/smallref.fa.* && rm -r temp_test && rm quatradis_out.profile
echo "ok"

echo -n "Checking 'tradis pipeline multiple' ... "
./scripts/tradis pipeline multiple --output_dir temp_test tests/data/mapper/fastq.txt tests/data/mapper/smallref.fa 2>&1 > /dev/null && rm -r temp_test work .nextflow && rm .nextflow.log*
echo "ok"
