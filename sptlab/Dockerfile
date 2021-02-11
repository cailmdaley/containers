ARG COSCINE_VERSION=latest
FROM cailmdaley/coscine:$COSCINE_VERSION as base

USER root
# dependencies for 3g-software
RUN apt-get update && apt-get install -y --no-install-recommends \
        libboost-all-dev \
        libfftw3-dev \
        libflac-dev \
        libgsl0-dev \
        libnetcdf-dev

USER $IMAGE_USER
RUN pip install --no-cache-dir \
        astropy \
        camb \
        healpy \
        h5py \
        numexpr
COPY --chown=1000 *.sh $SRC/bin/
