# This container will install Bio-Tradis from master
#
FROM python:3.10-slim-bullseye

# Install the dependancies
RUN apt-get update -qq && \
    apt-get install -y sudo bzip2 locales unzip wget && \
    apt-get install -y bwa minimap2 samtools smalt tabix

# Set locales (required for running in Singularity)
RUN   sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
      touch /usr/share/locale/locale.alias && \
      locale-gen
ENV   LANG     en_GB.UTF-8
ENV   LANGUAGE en_GB:en
ENV   LC_ALL   en_GB.UTF-8

# Add source code
ADD . biotradis2

# Set environment
ENV PATH /biotradis2/bin:$PATH
WORKDIR /work
