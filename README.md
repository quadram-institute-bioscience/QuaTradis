# QuaTradis (Quadram TraDis)

A set of tools to analyse the output from TraDIS analyses  

 
[![Build Status](https://img.shields.io/travis/com/quadram-institute-bioscience/QuaTradis/master)](https://app.travis-ci.com/github/quadram-institute-bioscience/QuaTradis) 
[![Docker Build Status](https://img.shields.io/docker/build/quadraminstitute/quatradis.svg)](https://hub.docker.com/r/quadraminstitute/quatradis) [![Docker Pulls](https://img.shields.io/docker/pulls/quadraminstitite/quatradis.svg)](https://hub.docker.com/r/quadraminstitute/quatradis)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/quatradis/README.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/QuaTradis/blob/master/LICENSE) [![status](https://img.shields.io/badge/Bioinformatics-10.1093-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btw022)  
<!--
[![codecov](https://codecov.io/gh/sanger-pathogens/bio-tradis/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/bio-tradis)
-->
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
- From bioconda, for a simple installation and user experience
- From docker, if you want to guarantee you are using the exact official version (although it comes with a little extra complexity to run)

Each installation method is described below.

### From source

[![Github last commit](https://img.shields.io/github/last-commit/quadram-institute-bioscience/quatradis)](https://github.com/quadram-institute-bioscience/QuaTradis)
[![Github last release](https://img.shields.io/github/release-date/quadram-institute-bioscience/quatradis)](https://github.com/quadram-institute-bioscience/QuaTradis)
[![Github downloads](https://img.shields.io/github/downloads/quadram-institute-bioscience/quatradis/total)](https://github.com/quadram-institute-bioscience/QuaTradis)
[![Github latest version](https://img.shields.io/github/tag/quadram-institute-bioscience/quatradis?sort=semver)](https://github.com/quadram-institute-bioscience/QuaTradis)  

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
  * python3-pip
  * r >= 3.6

```bash
sudo apt-get update -qq && \
sudo apt-get install -y bzip2 default-jre gcc locales make python3-pip r-base unzip wget && \
sudo apt-get install -y bwa minimap2 smalt
```
  
With the core packages installed we also need some python packages.  These can be installed using pip as described here:

```bash
pip3 install Bio cython pysam numpy pytest-cov semantic-version snakeviz
```

Finally, we also need a few R packages:

```bash
sudo Rscript -e "install.packages('BiocManager')" -e "BiocManager::install()" -e "BiocManager::install(c('edgeR','getopt', 'MASS'))"
```


Next move to a location to which you want to put the source code then clone the github repo and create a development build:

```buildoutcfg
git clone https://github.com/quadram-institute-bioscience/QuaTradis.git
cd QuaTradis
make dev
```
QuaTradis scripts should now be available.


#### Troubleshooting installations on old distributions

Older distributions of linux such as ubuntu 18 come packaged with R 3.4.  In this case R packages such as MASS and EdgeR will not install correctly with the commands mentioned above.  In this case we would recommend updating to R 3.6+.  In this snippet we show how to update R to 4.  However be cautious as this snippet will uninstall the current version of R, so be sure this is what you want to do, or figure out another way of running the versions side by side.

```bash
sudo apt install libssl-dev libcurl4-openssl-dev libxml2-dev
sudo apt remove r-base* --purge
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
sudo apt update
sudo apt install r-base
```


### Bioconda

[![Anaconda-Server Badge](https://anaconda.org/bioconda/quatradis/badges/version.svg)](https://anaconda.org/bioconda/quatradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/quatradis/badges/latest_release_date.svg)](https://anaconda.org/bioconda/quatradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/quatradis/badges/platforms.svg)](https://anaconda.org/bioconda/quatradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/quatradis/badges/downloads.svg)](https://anaconda.org/bioconda/quatradis)

Before installing quatradis via conda, first install conda and enable the bioconda channel as described in steps 1 and 2 
[here](http://bioconda.github.io/user/install.html#set-up-channels).

```
conda install -c bioconda quatradis
```

Note: If this is the first time you have installed a conda package on your system this may take some time.

### From Docker

[![Docker Build Status](https://img.shields.io/docker/build/quadraminstitute/quatradis.svg)](https://hub.docker.com/r/quadraminstitute/quatradis) 
[![Docker Pulls](https://img.shields.io/docker/pulls/quadraminstitite/quatradis.svg)](https://hub.docker.com/r/quadraminstitute/quatradis)
[![Docker Version](https://img.shields.io/docker/v/quadraminstitite/quatradis?sort=semver)](https://hub.docker.com/r/quadraminstitute/quatradis)  
[![Docker Image Size](https://img.shields.io/docker/image-size/quadraminstitite/quatradis?sort=semver)](https://hub.docker.com/r/quadraminstitute/quatradis)  

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

Executable scripts to carry out most of the listed functions are available in the `scripts`, although these should get
installed onto your path after installation.  All tradis functions (except those delivered by R scripts) have been condensed 
into a single executable named `tradis`.  Under this there are multiple sub-functions:

* `tags`
  * `add` - Generates a BAM file with tags added to read strings.
  * `check`- Prints 1 if tags are present in alignment file, prints 0 if not.
  * `filter` - Create a fastq file containing reads that match the supplied tag.
  * `remove` - Creates a fastq file containing reads with the supplied tag removed from the sequences.
* `plot`
  * `create` - Creates an gzipped insertion site plot from a set of alignments.
  * `combine` - Combines several plots into 1.
  * `analyse` - Takes genome annotation in embl format along with plot files produced by bacteria_tradis and generates tab-delimited files containing gene-wise annotations of insert sites and read counts.
* `pipeline`
  * `single` - Runs complete analysis, starting with a fastq file and produces mapped BAM files and plot files for each file in the given file list and a statistical summary of all files. Note that the -f option expects a text file containing a list of fastq files, one per line. This script can be run with or without supplying tags. 
  * `nextflow` - Same as single, put can process multiple fastq files in parallel using nextflow.  This is capable of distributing work over a cluster.
* `tradis_essentiality.R` - Takes a single tab-delimited file from tradis_gene_insert_sites to produce calls of gene essentiality. Also produces a number of diagnostic plots.
* `tradis_comparison.R` - Takes tab files to compare two growth conditions using edgeR. This analysis requires experimental replicates.

Note that default parameters are for comparative experiments, and will need to be modified for gene essentiality studies.

A help menu for each script can be accessed by running the script by adding with "--help".

### Docker

To use QuaTradis from docker then the command line gets a little more complicated.  Something like the command below,
should work however.  And once you find something that works for you, then a simple wrapper script can simplify execution
sigificantly.

```bash
    docker run --rm -it -u $(id -u ${USER}):$(id -g ${USER}) -v /home/ubuntu/data:/data quadraminstitute/quatradis <program and args>
```

To explain what is happening here you should familiarise yourself with [docker](https://docs.docker.com/engine/reference/commandline/run/) 
if you haven't already.  But a quick summary is:
- `--rm` Remove the container after it has completed.
- `-it` Run in an interactive terminal session
- `-u $(id -u ${USER}):$(id -g ${USER})` Run container as the current host user.  This stops all generated files being owned by root.
- `-v /home/ubuntu/data:/data` Mount the local directory `/home/ubuntu/data` to `/data` in the container.  Change the first half of the parameter as necessary to mount in your data files.

Replace `<program and args>` with whatever script you wish to use in this package, and provide whatever arguments make sense for that script.


## Contributing

Quatradis uses git for version control, with its public repo sited on [github](https://github.com/quadram-institute-bioscience/QuaTradis/tree/master).  
Contributions go through a standard feature branch workflow, so for most users this means that changes should be created 
in a feature branch off of the default master branch.  Changes are merged into `master` by submitting a 
[github pull request](https://github.com/quadram-institute-bioscience/QuaTradis/pulls).  

Before submitting the PR please check that your code changes pass all unit tests (see below) to save time.  Also ideally 
write some unit tests to cover your new functionality.  PRs require 1 approval from an admin as well as a successful travis 
build before the PR can be merged.

### Building quatradis from source

Quatradis is mostly python code, managed with [setuptools](https://pythonhosted.org/an_example_pypi_project/setuptools.html), 
so all the usual methods for building python via this method apply.  However, to make things simpler a 
[Makefile](Makefile) is provided which various targets for some common actions a developer might want to take.  For example 
if you want to build a development release type `make dev`.  If you would like to install into your local python environment
type `make install`.  A local docker image can be built by typing `make docker-build`.

### Running unit tests

The suite of unit tests can be run with `make test` from the top-level directory.  Alternatively you can run `pytest` 
directly from the `tests` directory.

### CI

[![Build Status](https://img.shields.io/travis/com/quadram-institute-bioscience/QuaTradis/master)](https://app.travis-ci.com/github/quadram-institute-bioscience/QuaTradis)

Continuous integration is delivered via [travis](https://app.travis-ci.com/github/quadram-institute-bioscience/QuaTradis).
The [travis pipeline](.travis.yml) is designed to build and test commits from PRs and the master branch.  For `master`
branch builds, travis will also publish the latest docker image to [dockerhub](https://hub.docker.com/r/quadraminstitute/quatradis) 
if that commit is tagged (see below). 

### Versioning

[![Github latest version](https://img.shields.io/github/tag/quadram-institute-bioscience/quatradis?sort=semver)](https://github.com/quadram-institute-bioscience/QuaTradis)  

When administrators (as defined by github) feel it is time to create a new release and bump the version there are 3 different
make targets to choose from depending on which kind of release is required:

- Major release: `make release_major`
- Minor release: `make release_minor`
- Patch release: `make release_patch`

For details of what each of these do see the [Makefile](Makefile) but this increases the dot release number
in the [VERSION](VERSION) file and commits this to the `master` branch and pushes it to github.  The new master commit 
is tagged appropriately with the new version.  This then triggers travis to build and publish a docker image to dockerhub.

If administrators need help figuring out which release number to bump then check the [semver2](https://semver.org/) guidelines. 


## License
QuaTradis is free software, licensed under [GPLv3](LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/quadram-institute-bioscience/QuaTradis/issues), or contact
the [developers](AUTHORS).

## Citation
If you use this software please cite:

"The TraDIS toolkit: sequencing and analysis for dense transposon mutant libraries", Barquist L, Mayho M, Cummins C, Cain AK, Boinett CJ, Page AJ, Langridge G, Quail MA, Keane JA, Parkhill J. Bioinformatics. 2016 Apr 1;32(7):1109-11. doi: [10.1093/bioinformatics/btw022](https://doi.org/10.1093/bioinformatics/btw022). Epub 2016 Jan 21.
