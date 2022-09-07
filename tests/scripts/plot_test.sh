#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir
DATA_DIR=tests/data/tisp

# Sanity test help messages
echo -n "Checking 'tradis plot' help messages ... "
./tradis plot --help > /dev/null
./tradis plot create --help  > /dev/null
./tradis plot combine --help > /dev/null
./tradis plot count --help > /dev/null
./tradis plot normalise --help > /dev/null
echo "ok"

echo -n "Checking 'tradis plot combine' ... "
./tradis plot combine $DATA_DIR/combine/zip_comb_list.txt > /dev/null 2>&1 && rm -r combined && rm zip_comb_list.stats
echo "ok"

echo -n "Checking 'tradis plot count' ... "
./tradis plot count -o temp_test $DATA_DIR/analyse/reference_BW25113_short.embl $DATA_DIR/analyse/controlLBrep?.insert_site_plot_short.gz > /dev/null 2>&1 && rm -r temp_test
echo "ok"

echo -n "Checking 'tradis plot create from_fastq' ... "
./tradis plot create from_fastq --output_dir temp_test --profile tests/data/mapper/test.fastq tests/data/mapper/smallref.fa > /dev/null 2>&1 && rm -f tests/data/mapper/smallref.fa.* && rm -r temp_test && rm tradis.profile
echo "ok"

echo -n "Checking 'tradis plot create from_alignments' ... "
./tradis plot create from_alignments --outfile temp_test/tradis.plot -m 20 $DATA_DIR/create/small_multi_sequence.bam > /dev/null 2>&1 && rm -r temp_test
echo "ok"

echo -n "Checking 'tradis plot normalise' ... "
./tradis plot normalise $DATA_DIR/normalise/fewinsertions $DATA_DIR/normalise/manyinsertions > /dev/null 2>&1 && rm $DATA_DIR/normalise/*.norm
echo "ok"