# This container will install QuaTradis from master
#
FROM debian:bullseye-slim

# Install the dependencies
RUN apt-get update -qq && \
    apt-get install -y sudo bzip2 curl default-jre gcc gzip locales make procps r-base unzip wget && \
    apt-get install -y bwa minimap2 smalt

# Install R dependencies
RUN Rscript -e "install.packages('BiocManager')" -e "BiocManager::install()" -e "BiocManager::install(c('edgeR','getopt', 'MASS'))"

# Install HTSlib to get bgzip installed as an executable
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && \
    make install

# Install snakemake
# Set locales (required for running in Singularity)
RUN sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
    touch /usr/share/locale/locale.alias && \
    locale-gen
ENV LANG     en_GB.UTF-8
ENV LANGUAGE en_GB:en
ENV LC_ALL   en_GB.UTF-8
ENV SHELL /bin/bash
ENV PATH /opt/conda/bin:${PATH}
RUN /bin/bash -c "curl -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh"
RUN /bin/bash -c "conda install -y -c conda-forge mamba && \
    mamba create -q -y -c conda-forge -c bioconda -n quatradis snakemake snakemake-minimal && \
    source activate quatradis && \
    mamba install -q -y -c conda-forge -c bioconda snakemake snakemake-minimal && \
    mamba install -q -y -c conda-forge singularity && \
    conda clean --all -y && \
    which python"
RUN echo "source activate quatradis" > ~/.bashrc
ENV PATH /opt/conda/envs/quatradis/bin:${PATH}
RUN pip install jinja2 pygments slacker oauth2client google-crc32c google-api-python-client google-cloud-storage

# Add source code
ADD . /quatradis
WORKDIR /quatradis

# Install python dependencies
RUN pip install -r requirements.txt
RUN pip install .[test]

# Install quatradis
RUN pip install .

# Set default working dir
WORKDIR /work
