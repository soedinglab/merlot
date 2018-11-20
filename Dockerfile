FROM frolvlad/alpine-python3

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

LABEL maintainer="npapado@mpibpc.mpg.de"

# install python dependencies
RUN pip install scipy
RUN pip install cython
RUN pip install git+https://github.com/soedinglab/csgraph_mod

# Use R 3.4 since not all dependencies are updated to 3.5
FROM rocker/tidyverse:3.4

# install bash dependencies
RUN apt-get update
RUN apt-get install -y libudunits2-dev
RUN apt-get install -y mesa-common-dev
RUN apt-get install -y libglu1-mesa-dev

# install R packages
RUN ["Rscript", "inst/scripts/install_R.R"]

# run the selector script that will fire up other scripts
ENTRYPOINT ["inst/scripts/selector.sh"]