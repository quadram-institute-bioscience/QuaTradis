#!/bin/bash

set -e

echo -n "Checking 'tradis' help message ... "
./tradis --help > /dev/null
echo "ok"

./tests/scripts/utils_test.sh
./tests/scripts/tags_test.sh
./tests/scripts/plot_test.sh
./tests/scripts/comparison_test.sh
./tests/scripts/pipeline_test.sh
./tests/scripts/R_test.sh