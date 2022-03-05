#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir

# Sanity test help messages
echo -n "Checking 'tradis utils' help messages ... "
./scripts/tradis utils --help > /dev/null
./scripts/tradis utils index --help > /dev/null
echo "ok"

echo -n "Checking 'tradis utils index' ... "
./scripts/tradis utils index -a bwa tests/data/mapper/smallref.fa testbwaref && rm tests/data/mapper/smallref.fa.*
./scripts/tradis utils index -a smalt tests/data/mapper/smallref.fa testsmaltref && rm testsmaltref.*
./scripts/tradis utils index -a minimap2 tests/data/mapper/smallref.fa testminimapref && rm testminimapref
echo "ok"
