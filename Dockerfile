FROM python:3.6

LABEL maintainer="npapado@mpibpc.mpg.de"

# install python dependencies
RUN pip install scipy
RUN pip install cython
RUN pip install git+https://github.com/soedinglab/csgraph_mod

# Use R 3.4 since not all dependencies are updated to 3.5
FROM rocker/tidyverse:3.5

# Set the working directory to /app
WORKDIR /data

# Copy the current directory contents into the container at /app
COPY . /app

# install bash dependencies
RUN apt-get update
RUN apt-get install -y libudunits2-dev
RUN apt-get install -y mesa-common-dev
RUN apt-get install -y libglu1-mesa-dev

# install R packages
RUN ["Rscript", "/app/inst/scripts/install_R.R"]
COPY inst/scripts/* /app/

# run the selector script that will fire up other scripts
ENTRYPOINT ["/app/docker-entrypoint.sh"]
CMD ["--help"]