# QuaTradis (Quadram TraDis)

A set of tools to analyse the output from TraDIS analyses  

<!--
[![Build Status](https://travis-ci.org/sanger-pathogens/Bio-Tradis.svg?branch=master)](https://travis-ci.org/sanger-pathogens/Bio-Tradis)  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/Bio-Tradis/blob/master/software_license)  
[![status](https://img.shields.io/badge/Bioinformatics-10.1093-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btw022)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/biotradis/README.html)   
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/biotradis)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/bio-tradis.svg)](https://hub.docker.com/r/sangerpathogens/bio-tradis)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/bio-tradis.svg)](https://hub.docker.com/r/sangerpathogens/bio-tradis)  
[![codecov](https://codecov.io/gh/sanger-pathogens/bio-tradis/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/bio-tradis)
-->

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [Bioconda](#bioconda)
    * [Docker](#docker)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
    * [Scripts](#scripts)
    * [Analysis Scripts](#analysis-scripts)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction 
The QuaTradis pipeline provides software utilities for the processing, mapping, and analysis of transposon insertion sequencing data. The pipeline was designed with the data from the TraDIS sequencing protocol in mind, but should work with a variety of transposon insertion sequencing protocols as long as they produce data in the expected format.

For more information on the TraDIS method, see http://bioinformatics.oxfordjournals.org/content/32/7/1109 and http://genome.cshlp.org/content/19/12/2308.

## Installation
QuaTradis has the following dependencies:

### Required dependencies
* bwa
* smalt
* samtools
* tabix

There are a number of ways to install QuaTradis and details are provided below. If you encounter an issue when installing QuaTradis please contact your local system administrator. 

### Bioconda
Install conda and enable the bioconda channel.
<!--
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/version.svg)](https://anaconda.org/bioconda/biotradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/latest_release_date.svg)](https://anaconda.org/bioconda/biotradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/platforms.svg)](https://anaconda.org/bioconda/biotradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/downloads.svg)](https://anaconda.org/bioconda/biotradis)
-->
```
conda install -c bioconda quatradis=xxx
```

### Docker
QuaTradis can be run in a Docker container. First install Docker, then pull the QuaTradis image from dockerhub:

    docker pull quadraminstitute/quatradis

To use QuaTradis use a command like this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:

    docker run --rm -it -v /home/ubuntu/data:/data quadraminstitute/quatradis bacteria_tradis -h


## Usage

QuaTradis provides functionality to:
* detect TraDIS tags in a BAM file
* add the tags to the reads
* filter reads in a FastQ file containing a user defined tag
* remove tags
* map to a reference genome
* create an insertion site plot file
  
The functions are available as standalone scripts or as perl modules.

### Scripts
Executable scripts to carry out most of the listed functions are available in the `bin`:

* `check_tradis_tags` - Prints 1 if tags are present in alignment file, prints 0 if not.
* `add_tradis_tags` - Generates a BAM file with tags added to read strings.
* `filter_tradis_tags` - Create a fastq file containing reads that match the supplied tag
* `remove_tradis_tags` - Creates a fastq file containing reads with the supplied tag removed from the sequences
* `tradis_plot` - Creates an gzipped insertion site plot
* `bacteria_tradis` - Runs complete analysis, starting with a fastq file and produces mapped BAM files and plot files for each file in the given file list and a statistical summary of all files. Note that the -f option expects a text file containing a list of fastq files, one per line. This script can be run with or without supplying tags. 
* `tradis_gene_insert_sites` - Takes genome annotation in embl format along with plot files produced by bacteria_tradis and generates tab-delimited files containing gene-wise annotations of insert sites and read counts.
* `tradis_essentiality.R` - Takes a single tab-delimited file from tradis_gene_insert_sites to produce calls of gene essentiality. Also produces a number of diagnostic plots.
* `tradis_comparison.R` - Takes tab files to compare two growth conditions using edgeR. This analysis requires experimental replicates.

Note that default parameters are for comparative experiments, and will need to be modified for gene essentiality studies.

A help menu for each script can be accessed by running the script by adding with "--help".

## Contributing

All changes to Quatradis should be made in a separate git branch from master, and then integrated into the master branch using
a [github pull request](https://github.com/quadram-institute-bioscience/QuaTradis/pulls).  Before submitting the PR
please check that your code changes pass all unit tests (see below).  Also ideally write some unit tests to cover your
new functionality.  PRs currently require 1 approval and a successful travis build.

### Running unit tests

The test can be run with `pytest` from the `tests` directory.  
Alternatively you can use the make target from the top-level
directory:  `make test`  

### Travis CI

Continuous integration is delivered via [travis](https://app.travis-ci.com/github/quadram-institute-bioscience/QuaTradis).
The [travis pipeline](.travis.yml) is designed to build and test all commit on all branches pushed to github.  For master
builds, which only take place after a successful PR, travis will also publish the latest docker image to 
[dockerhub](https://hub.docker.com/r/quadraminstitute/quatradis). 

### Versioning

When administrators of the 


## License
QuaTradis is free software, licensed under [GPLv3](LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/quadram-institute-bioscience/QuaTradis/issues).

## Citation
If you use this software please cite:

"The TraDIS toolkit: sequencing and analysis for dense transposon mutant libraries", Barquist L, Mayho M, Cummins C, Cain AK, Boinett CJ, Page AJ, Langridge G, Quail MA, Keane JA, Parkhill J. Bioinformatics. 2016 Apr 1;32(7):1109-11. doi: [10.1093/bioinformatics/btw022](https://doi.org/10.1093/bioinformatics/btw022). Epub 2016 Jan 21.
