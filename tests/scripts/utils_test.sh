#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir
DATA_DIR=tests/data
MAPPER_DATA_DIR=$DATA_DIR/util/mapper
EMBL_DATA_DIR=$DATA_DIR/embl
ESS_DATA_DIR=$DATA_DIR/essentiality/input

# Sanity test help messages
echo -n "Checking 'tradis utils' help messages ... "
./tradis utils --help > /dev/null
./tradis utils index --help > /dev/null
./tradis utils annotation --help > /dev/null
./tradis utils artemis_project --help > /dev/null
./tradis utils gene_reports --help > /dev/null
./tradis utils prepare_embl --help > /dev/null
./tradis utils essentiality --help > /dev/null
echo "ok"

echo -n "Checking 'tradis utils index' ... "
./tradis utils index -a bwa $MAPPER_DATA_DIR/smallref.fa testbwaref && rm $MAPPER_DATA_DIR/smallref.fa.*
./tradis utils index -a smalt $MAPPER_DATA_DIR/smallref.fa testsmaltref && rm testsmaltref.*
./tradis utils index -a minimap2 $MAPPER_DATA_DIR/smallref.fa testminimapref && rm testminimapref
echo "ok"

echo -n "Checking 'tradis utils prepare_embl' ... "
./tradis utils prepare_embl --output valid_prepared.embl --minimum_threshold=1 --window_size=4 --window_interval=2 \
  --prime_feature_size=100 $EMBL_DATA_DIR/prepareinputfiles/valid && rm valid_prepared.embl
./tradis utils prepare_embl --output valid_prepared.embl --minimum_threshold=1 --window_size=4 --window_interval=2 \
  --prime_feature_size=100 --emblfile=$EMBL_DATA_DIR/expandgenes/one_gene $EMBL_DATA_DIR/prepareinputfiles/valid && rm valid_prepared.embl
echo "ok"

echo -n "Checking 'tradis utils essentiality' ... "
./tradis utils essentiality --output_dir=output_essential --minimum_threshold=1 $ESS_DATA_DIR/small_case.insert_site_plot.gz $ESS_DATA_DIR/annotation.embl && rm -r output_essential
echo "ok"
