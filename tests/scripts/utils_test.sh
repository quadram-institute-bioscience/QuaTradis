#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir
DATA_DIR=tests/data
MAPPER_DATA_DIR=$DATA_DIR/util/mapper

# Sanity test help messages
echo -n "Checking 'tradis utils' help messages ... "
./tradis utils --help > /dev/null
./tradis utils index --help > /dev/null
./tradis utils annotation --help > /dev/null
./tradis utils artemis_project --help > /dev/null
./tradis utils extract_names --help > /dev/null
./tradis utils tags --help > /dev/null
echo "ok"

echo -n "Checking 'tradis utils index' ... "
./tradis utils index -a bwa $MAPPER_DATA_DIR/smallref.fa testbwaref
rm $MAPPER_DATA_DIR/smallref.fa.*
./tradis utils index -a smalt $MAPPER_DATA_DIR/smallref.fa testsmaltref
rm testsmaltref.*
./tradis utils index -a minimap2 $MAPPER_DATA_DIR/smallref.fa testminimapref
rm testminimapref
echo "ok"


