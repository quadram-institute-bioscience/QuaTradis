#!/bin/bash

# Exit immediately on error
#set -e

# Sanity test help messages
echo -n "Checking 'tradis_comparison.R script' help message ... "
./quatradis/essentiality/tradis_comparison.R --help > /dev/null 2>&1
echo "ok"

echo -n "Checking 'tradis_essentiality.R script' help message ... "
./quatradis/essentiality/tradis_essentiality.R > /dev/null 2>&1
echo "ok"