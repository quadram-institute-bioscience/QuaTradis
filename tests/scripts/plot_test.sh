#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir

# Sanity test help messages
echo -n "Checking 'tradis plot' help messages ... "
./scripts/tradis plot --help > /dev/null
./scripts/tradis plot create --help  > /dev/null
./scripts/tradis plot combine --help > /dev/null
./scripts/tradis plot analyse --help > /dev/null
echo "ok"

echo -n "Checking 'tradis plot create' ... "
./scripts/tradis plot create --outfile temp_test/tradis.plot -m 20 tests/data/isp_create/small_multi_sequence.bam 2>&1 && rm -r temp_test
echo "ok"

echo -n "Checking 'tradis plot combine' ... "
./scripts/tradis plot combine tests/data/isp_combine/zip_comb_list.txt 2>&1 && rm -r combined && rm zip_comb_list.stats
echo "ok"

echo -n "Checking 'tradis plot analyse' ... "
./scripts/tradis plot analyse -o temp_test tests/data/isp_analyse/reference_BW25113_short.embl tests/data/isp_analyse/controlLBrep?.insert_site_plot_short.gz && rm -r temp_test
echo "ok"
