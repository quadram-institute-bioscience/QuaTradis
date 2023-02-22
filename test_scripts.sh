#!/bin/bash
./tradis --help > /dev/null
./tests/scripts/tags_test.sh
./tests/scripts/plot_test.sh
./tests/scripts/comparison_test.sh
./tests/scripts/utils_test.sh
./tests/scripts/pipeline_test.sh
./tests/scripts/R_test.sh