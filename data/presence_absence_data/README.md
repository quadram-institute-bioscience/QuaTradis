# Example data for tradis ca presence_absence function

This folder contains some example data to run the presence_absence workflow. The data is from the Triclosan dataset presented in Yasir et al. 2019 (https://www.biorxiv.org/content/10.1101/596833v1). 

## Input files
List of files: 
* reference_BW25113.embl - E-coli BW25113 Genome annotation 
* gene_report_0008mgL.csv - gene report for 0.008mg Triclosan concentration
* gene_report_0015mgL.csv - gene report for 0.015mg Triclosan concentration
* gene_report_003mgL.csv - gene report for 0.03mg Triclosan concentration
* gene_report_006mgL.csv - gene report for 0.06mg Triclosan concentration
* gene_report_0125mgL.csv - gene report for 0.125mg Triclosan concentration
* gene_report_025mgL.csv - gene report for 0.25mg Triclosan concentration
* gene_report_05mgL.csv - gene report for 0.5mg Triclosan concentration
* gene_report_1mgL.csv - gene report for 1mg Triclosan concentration


## Commands

This command can be run from this working directory, without any changes, assuming QuaTraDIS has been successfully installed. It is a good way to see what the output looks like and to test if the installation has been performed correctly.
```
tradis ca presence_absence reference_BW25113.embl gene_report_0008mgL.csv gene_report_0015mgL.csv gene_report_003mgL.csv gene_report_006mgL.csv gene_report_0125mgL.csv gene_report_025mgL.csv gene_report_05mgL.csv gene_report_1mgL.csv
```
This command only takes a few seconds to run on a standard laptop.

Alternatively if you don't have QuaTraDIS installed, you can copy and paste the following docker command (must have docker installed), which will download an installation of QuaTraDIS and run it without any changes.
```
#TODO test this
docker run --rm -it -v $PWD/tradis_output:/work quadraminstitute/quatradis:latest tradis ca presence_absence QuaTradis/data/presence_absence_data/reference_BW25113.embl QuaTradis/data/presence_absence_data/gene_report_0008mgL.csv QuaTradis/data/presence_absence_data/gene_report_0015mgL.csv QuaTradis/data/presence_absence_data/gene_report_003mgL.csv QuaTradis/data/presence_absence_data/gene_report_006mgL.csv QuaTradis/data/presence_absence_data/gene_report_0125mgL.csv QuaTradis/data/presence_absence_data/gene_report_025mgL.csv QuaTradis/data/presence_absence_data/gene_report_05mgL.csv QuaTradis/data/presence_absence_data/gene_report_1mgL.csv
```