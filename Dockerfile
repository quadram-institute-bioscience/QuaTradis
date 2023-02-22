# This container will install QuaTradis from master
#
FROM sbastkowski/quatradis-base:latest

# Install python dependencies
ADD requirements.txt /home/quatradis/src/requirements.txt
WORKDIR /home/quatradis/src
RUN pip install -r requirements.txt

# Add full source code
ADD . /home/quatradis/src

# Install quatradis
RUN pip install .

# Run unit tests
RUN cd tests && pytest --doctest-modules --disable-warnings

# Run script test
RUN ./test_scripts.sh

# Set default working dir
WORKDIR /work
