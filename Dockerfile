# This container will install QuaTradis from master
#
FROM python:3.10-slim-bullseye

# Install the dependencies
RUN apt-get update -qq && \
    apt-get install -y sudo bzip2 gcc default-jre locales unzip wget make && \
    apt-get install -y bwa minimap2 smalt

# Install nextflow
RUN wget -qO- https://get.nextflow.io | bash && \
	chmod +x nextflow && \
	mv nextflow /usr/local/bin

RUN pip install Bio pysam numpy pytest semantic-version snakeviz

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
RUN python3 setup.py install

# Set environment
ENV PATH /quatradis/bin:$PATH
WORKDIR /work
