#!/bin/bash

# Exit immediately on error
set -e

# This script assumes it is run from project root dir
DATA_DIR=tests/data
PIPELINE_DATA_DIR=tests/data/pipeline/input

# Sanity test help messages
echo -n "Checking 'tradis compare' help messages ... "
./tradis compare --help > /dev/null
./tradis compare figures --help  > /dev/null
./tradis compare presence_absence --help  > /dev/null
./tradis compare insertion_sites --help  > /dev/null
./tradis compare gene_report --help  > /dev/null
./tradis compare split --help  > /dev/null
./tradis compare essentiality --help  > /dev/null
./tradis compare essentiality_analysis --help  > /dev/null
./tradis compare prepare_embl --help  > /dev/null
echo "ok"

echo -n "Checking 'tradis compare presence_abscence' ... "
PA_DATA=$DATA_DIR/comparison/presenceabsence
./tradis compare presence_absence -o $PA_DATA/output \
  $PA_DATA/reference.embl $PA_DATA/ctrl.csv $PA_DATA/gent1.csv $PA_DATA/gent5.csv $PA_DATA/gent25.csv> /dev/null 2>&1
cmp $PA_DATA/expected_logfc.csv $PA_DATA/output/all_logfc.csv
rm -r $PA_DATA/output
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

echo -n "Checking 'tradis compare essentiality' ... "
ESS_DATA_DIR=$DATA_DIR/comparison/essentiality
./tradis compare essentiality $ESS_DATA_DIR/combined.count.tsv
ESS_OUT=$ESS_DATA_DIR/combined.count.tsv.essen.csv
if [ ! -f "$ESS_OUT" ]; then
    echo "tradis compare essentiality test failed.  $ESS_OUT does not exist."
    exit 1
fi
rm $ESS_DATA_DIR/combined.count.tsv.*
echo "ok"

echo -n "Checking 'tradis compare prepare_embl' under comparison_test.sh ... "
EMBL_DATA_DIR=$DATA_DIR/embl
./tradis compare prepare_embl --output valid_prepared.embl --minimum_threshold=1 --window_size=4 --window_interval=2 \
  --plotfile $PIPELINE_DATA_DIR/small_case.insert_site_plot.gz $PIPELINE_DATA_DIR/small_case_2.insert_site_plot.gz $PIPELINE_DATA_DIR/small_control.insert_site_plot.gz $PIPELINE_DATA_DIR/small_control_high_insertions.insert_site_plot.gz --prime_feature_size=100 --emblfile=$EMBL_DATA_DIR/expandgenes/one_gene --dynamic_window
echo "ok prepare embl dynamic window tested successfully"
rm valid_prepared.embl
./tradis compare prepare_embl --output valid_prepared.embl --minimum_threshold=1 --window_size=4 --window_interval=2 \
  --plotfile $PIPELINE_DATA_DIR/small_case.insert_site_plot.gz $PIPELINE_DATA_DIR/small_case_2.insert_site_plot.gz $PIPELINE_DATA_DIR/small_control.insert_site_plot.gz $PIPELINE_DATA_DIR/small_control_high_insertions.insert_site_plot.gz --prime_feature_size=100 --emblfile=$EMBL_DATA_DIR/expandgenes/one_gene
echo "ok prepare embl without dynamic_window tested successfully"
rm valid_prepared.embl
echo "ok-final"