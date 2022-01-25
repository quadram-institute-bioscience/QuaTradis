# QuaTradis (Quadram TraDis)

A set of tools to analyse the output from TraDIS analyses  

 
[![Build Status](https://img.shields.io/travis/com/quadram-institute-bioscience/QuaTradis/master)](https://app.travis-ci.com/github/quadram-institute-bioscience/QuaTradis) [![Docker Build Status](https://img.shields.io/docker/build/quadraminstitute/quatradis.svg)](https://hub.docker.com/r/quadraminstitute/quatradis) [![Docker Pulls](https://img.shields.io/docker/pulls/quadraminstitite/quatradis.svg)](https://hub.docker.com/r/quadraminstitute/quatradis)  
<!--
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/quatradis/README.html)
-->
<!--
[![codecov](https://codecov.io/gh/sanger-pathogens/bio-tradis/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/bio-tradis)
-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/QuaTradis/blob/master/LICENSE) [![status](https://img.shields.io/badge/Bioinformatics-10.1093-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btw022)  

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [From source](#from-source)
    <!--* [From Bioconda](#from-bioconda)-->
    * [From docker](#from-docker)
  * [Usage](#usage)
    * [Scripts](#scripts)
    * [Docker](#docker)
  * [Contributing](#contributing)
    * [Building QuaTradis](#building-quatradis-from-source)
    * [Running unit tests](#running-unit-tests)
    * [Versioning](#versioning)
    * [CI](#ci)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction 

The QuaTradis pipeline provides software utilities for the processing, mapping, and analysis of transposon insertion 
sequencing data. The pipeline was designed with the data from the TraDIS sequencing protocol in mind, but should work 
with a variety of transposon insertion sequencing protocols as long as they produce data in the expected format.

QuaTradis provides functionality to:
* detect TraDIS tags in a BAM file
* add the tags to the reads
* filter reads in a FastQ file containing a user defined tag
* remove tags
* map to a reference genome
* create an insertion site plot file

For more information on the TraDIS method, see http://bioinformatics.oxfordjournals.org/content/32/7/1109 and 
http://genome.cshlp.org/content/19/12/2308.

## Installation

QuaTradis is designed to run on debian linux distros.  It may also run on other linux distros or MacOS but this has
not been tested.

There are multiple paths to installing QuaTradis depending on your use case.

- From source if you want to extend or debug functions
<!--- From bioconda, for a simple installation and user experience-->
- From docker, if you want to guarentee you are using the exact official version (although comes with a little extra complexity to run)

Each installation method is described below.

### From source

QuaTradis has the following dependencies so install these first.  We provide some guidelines on how to install on a debian
based system here, but if you are using another OS the actual commands may vary.  Also you may wish to install some of
these packages through alternate methods, e.g. from source.  In which case you'll need to lookup and follow the 
instructions for each package instead of running the commands below.

* System packages:
  * git
  * bwa
  * smalt
  * make
  * gcc
  * python3

```bash
apt-get update -qq && \
apt-get install -y sudo bzip2 gcc locales make python3 unzip wget && \
apt-get install -y bwa minimap2 smalt
```
  
With the core packages installed we also need some python packages.  These can be installed using pip as described here:

* Python packages:
  * Bio
  * pysam
  * numpy
  * semantic-version
  * snakeviz

```bash
pip install Bio pysam numpy pytest semantic-version snakeviz
```

Next move to a location to which you want to put the source code then clone the github repo and create a development build:

```buildoutcfg
git clone https://github.com/quadram-institute-bioscience/QuaTradis.git
cd QuaTradis
make dev
```
QuaTradis scripts should now be available.

<!--
### Bioconda
Install conda and enable the bioconda channel.

[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/version.svg)](https://anaconda.org/bioconda/biotradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/latest_release_date.svg)](https://anaconda.org/bioconda/biotradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/platforms.svg)](https://anaconda.org/bioconda/biotradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/biotradis/badges/downloads.svg)](https://anaconda.org/bioconda/biotradis)

```
conda install -c bioconda quatradis=xxx
```
-->
### From Docker
QuaTradis can be run in a Docker container. First install Docker, then pull the QuaTradis image from dockerhub:

    docker pull quadraminstitute/quatradis

if you wish to use a specific version, check what tags are available in [dockerhub](https://hub.docker.com/r/quadraminstitute/quatradis) 
and add `:<tagversion>` to the previous command (replacing with whatever version you wish to use).


## Usage

QuaTradis has a number of scripts which can be used from a shell to execute the functionality they wish.  More details 
of each script are described below.

In addition, the functions behind these scripts are available for use in the quatradis python package.  Details of using 
these functions is beyond this document, and left to individual developers to figure out.  Please contact the core
QuaTradis developers 

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

### Docker

To use QuaTradis from docker then the command line gets a little more complicated.  Something like the command below,
should work however.  And once you find something that works for you, then a simple wrapper script can simplify execution
sigificantly.

```bash
    docker run --rm -it -u $(id -u ${USER}):$(id -g ${USER}) -v /home/ubuntu/data:/data quadraminstitute/quatradis bacteria_tradis <program args>
```

To explain what is happening here you should familiarise yourself with [docker](https://docs.docker.com/engine/reference/commandline/run/) 
if you haven't already.  But a quick summary is:
- `--rm` Remove the container after it has completed.
- `-it` Run in an interactive terminal session
- `-u $(id -u ${USER}):$(id -g ${USER})` Run container as the current host user.  This stops all generated files being owned by root.
- `-v /home/ubuntu/data:/data` Mount the local directory `/home/ubuntu/data` to `/data` in the container.  Change the first half of the parameter as necessary to mount in your data files.

Feel free to replace the `bacteria_tradis` script with whatever script you wish to use in this package, and provide whatever
argument make sense for that script.


## Contributing

Quatradis uses git for version control, with its public repo sited on [github](https://github.com/quadram-institute-bioscience/QuaTradis/tree/master).  
Quatradis uses the [gitflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow)
branching workflow.  The reason for this is it allows contributions to go through a standard feature branch workflow, while
allows Quatradis admins to be selective about at which point in the development trunk branch a versioned release is made.
For most users this means that changes should be created in a feature branch off of the develop branch.  Changes should 
be submitted via [github pull request](https://github.com/quadram-institute-bioscience/QuaTradis/pulls) to the `develop` branch.  

Before submitting the PR please check that your code changes pass all unit tests (see below) to save time.  Also ideally 
write some unit tests to cover your new functionality.  PRs currently require 1 approval and a successful travis build.

### Building quatradis from source

Quatradis is mostly python code, managed with [setuptools](https://pythonhosted.org/an_example_pypi_project/setuptools.html), 
so all the usual methods for building python via this method apply.  However, to make things simpler a 
[Makefile](Makefile) is provided which various targets for some common actions a developer might want to take.  For example 
if you want to build a development release type `make dev`.  If you would like to install into your local python environment
type `make install`.  A local docker image can be built by typing `make docker-build`.

### Running unit tests

The suite of unit tests can be run with `make test` from the top-level directory.  Alternatively you can run `pytest` 
directly from the `tests` directory.

### Versioning

When administrators (as defined by github) feel it is time to create a new release and bump the version there are 3 different
make targets to choose from depending on which kind of release is required:

- Major release: `make release_major`
- Minor release: `make release_minor`
- Patch release: `make release_patch`

For details of what each of these do see the [Makefile](Makefile) but this increases the dot release number
in the [VERSION](VERSION) file and commits this to the `develop` branch.  It then merges this commit into `master` and pushes both 
branches to github.  The new master commit is also tagged with the new version.

If administrators need help figuring out which release to make check the [semver2](https://semver.org/) guidelines. 

### CI

Continuous integration is delivered via [travis](https://app.travis-ci.com/github/quadram-institute-bioscience/QuaTradis).
The [travis pipeline](.travis.yml) is designed to build and test all commits on all branches pushed to github.  For master
builds (i.e. new releases), travis will also publish the latest docker image to 
[dockerhub](https://hub.docker.com/r/quadraminstitute/quatradis). 


## License
QuaTradis is free software, licensed under [GPLv3](LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/quadram-institute-bioscience/QuaTradis/issues), or contact
the [developers](AUTHORS).

## Citation
If you use this software please cite:

"The TraDIS toolkit: sequencing and analysis for dense transposon mutant libraries", Barquist L, Mayho M, Cummins C, Cain AK, Boinett CJ, Page AJ, Langridge G, Quail MA, Keane JA, Parkhill J. Bioinformatics. 2016 Apr 1;32(7):1109-11. doi: [10.1093/bioinformatics/btw022](https://doi.org/10.1093/bioinformatics/btw022). Epub 2016 Jan 21.
