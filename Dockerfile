FROM alpine:3.8

LABEL maintainer="npapado@mpibpc.mpg.de"

WORKDIR /code

RUN apk add --no-cache --update\
    python3-dev\
    g++\
    bash\
    cmake\
    make\
    git\
    openblas-dev

RUN pip3 install --no-cache-dir numpy scipy cython pandas sklearn
RUN pip3 install --no-cache-dir git+https://github.com/soedinglab/csgraph_mod.git

# # make the directory that communicates with the outside world
WORKDIR /data

COPY inst/python/ScaffoldTree.py /code
# run the selector script that will fire up other scripts
ENTRYPOINT ["python3", "/code/ScaffoldTree.py"]