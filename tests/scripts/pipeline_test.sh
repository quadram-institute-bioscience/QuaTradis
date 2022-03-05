#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir

# Sanity test help messages
echo -n "Checking 'tradis pipeline' help messages ... "
./scripts/tradis pipeline --help | grep pipeline > /dev/null
./scripts/tradis pipeline single --help | grep pipeline > /dev/null
./scripts/tradis pipeline multiple --help | grep pipeline > /dev/null
echo "ok"

echo -n "Checking 'tradis pipeline single' ... "
./scripts/tradis pipeline single --output_dir temp_test tests/data/mapper/test.fastq tests/data/mapper/smallref.fa && rm tests/data/mapper/smallref.fa.* && rm -r temp_test
echo "ok"

echo -n "Checking 'tradis pipeline multiple' ... "
./scripts/tradis pipeline multiple --output_dir temp_test tests/data/mapper/fastq.txt tests/data/mapper/smallref.fa 2>&1 > /dev/null && rm -r temp_test work .nextflow && rm .nextflow.log*
echo "ok"
