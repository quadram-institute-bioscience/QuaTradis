#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir
ESSENTIALITY_DATA_DIR=tests/data/essentiality

# Sanity test help messages
echo -n "Checking 'tradis pipeline' help messages ... "
./tradis pipeline --help > /dev/null
./tradis pipeline create_plots --help > /dev/null
./tradis pipeline compare --help > /dev/null
echo "ok"

echo -n "Checking 'tradis pipeline create_plots' ... "
./tradis pipeline create_plots --output_dir temp_test tests/data/util/mapper/fastq.txt tests/data/util/mapper/smallref.fa > /dev/null 2>&1 && rm -f tests/data/util/mapper/smallref.fa.* && rm -r temp_test .snakemake
echo "ok"

echo -n "Checking 'tradis pipeline compare' ... "
./tradis pipeline compare --output_dir temp_test ESSENTIALITY_DATA_DIR/small_case.insert_site_plot.gz ESSENTIALITY_DATA_DIR/small_control.insert_site_plot.gz > /dev/null 2>&1 && rm -r temp_test .snakemake
echo "ok"