#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir

# Sanity test help messages
echo -n "Checking 'tradis compare' help messages ... "
./tradis compare --help > /dev/null
./tradis compare figures --help  > /dev/null
./tradis compare presence_absence --help  > /dev/null
./tradis compare analyse --help  > /dev/null
./tradis compare logfc_plot --help  > /dev/null
./tradis compare gene_report --help  > /dev/null
echo "ok"

echo -n "Checking 'tradis compare presence_abscence' ... "
DATA=data/presence_absence_data && \
./tradis compare presence_absence \
  $DATA/reference_BW25113.embl \
  $DATA/gene_report_0008mgL.csv $DATA/gene_report_0015mgL.csv \
  $DATA/gene_report_003mgL.csv $DATA/gene_report_006mgL.csv \
  $DATA/gene_report_0125mgL.csv $DATA/gene_report_025mgL.csv \
  $DATA/gene_report_05mgL.csv $DATA/gene_report_1mgL.csv > /dev/null 2>&1 && \
  rm -r output
echo "ok"

echo -n "Checking 'tradis compare analyse' ... "
DATA=data/comparative_analysis_data && \
./tradis compare analyse -v -a \
  $DATA/reference_BW25113_short.embl \
  $DATA/025mgLTricRep1.insert_site_plot_short.gz  \
  $DATA/controlLBrep1.insert_site_plot_short.gz > /dev/null 2>&1 && \
  rm -r output
echo "ok"

#echo -n "Checking 'tradis compare logfc_plot' ... "
#DATA=data/comparative_analysis_data && \
#./tradis compare logfc_plot -v -a \
#  $DATA/reference_BW25113_short.embl \
#  $DATA/025mgLTricRep1.insert_site_plot_short.gz  \
#  $DATA/controlLBrep1.insert_site_plot_short.gz > /dev/null 2>&1 && \
#  rm -r output
#echo "ok"

#echo -n "Checking 'tradis compare gene_report' ... "
#DATA=tests/data/comparison && \
#./tradis compare gene_report -v \
#  $DATA/analyse $DATA/annotation.embl > /dev/null 2>&1 && \
#  rm -r output
#echo "ok"