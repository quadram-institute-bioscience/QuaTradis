# This container will install QuaTradis from master
#
FROM python:3.10-slim-bullseye

# Install the dependencies
RUN apt-get update -qq && \
    apt-get install -y sudo bzip2 default-jre gcc gzip locales make procps r-base unzip wget && \
    apt-get install -y bwa minimap2 smalt

# Install nextflow
RUN wget -qO- https://get.nextflow.io | bash && \
	chmod +x nextflow && \
	mv nextflow /usr/local/bin

# Install HTSlib to get bgzip installed as an executable
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && \
    make install

# Install R dependencies
RUN Rscript -e "install.packages('BiocManager')" -e "BiocManager::install()" -e "BiocManager::install(c('edgeR','getopt', 'MASS'))"

# Set locales (required for running in Singularity)
RUN   sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
      touch /usr/share/locale/locale.alias && \
      locale-gen
ENV   LANG     en_GB.UTF-8
ENV   LANGUAGE en_GB:en
ENV   LC_ALL   en_GB.UTF-8

# Add source code
ADD . quatradis
WORKDIR /quatradis

RUN pip install -r requirements.txt
RUN python3 setup.py install

# Set environment
ENV PATH /quatradis/bin:$PATH
WORKDIR /work
