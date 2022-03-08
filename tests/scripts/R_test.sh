#!/bin/bash

# Exit immediately on error
#set -e

# Sanity test help messages
echo -n "Checking 'tradis_comparison.R script' help message ... "
./scripts/tradis_comparison.R --help 2>&1 > /dev/null
echo "ok"

echo -n "Checking 'tradis_essentiality.R script' help message ... "
./scripts/tradis_essentiality.R 2>&1 > /dev/null
echo "ok"