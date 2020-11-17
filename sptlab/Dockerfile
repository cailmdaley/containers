ARG COSCINE_VERSION=latest
FROM cailmdaley/coscine-tweaked:$COSCINE_VERSION

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        git-lfs \
        libboost-all-dev \
        libfftw3-dev \
        libflac-dev \
        libgsl0-dev \
        libnetcdf-dev \
        openssh-server

USER $IMAGE_USER

RUN pip install --no-cache-dir \
        astropy \
        camb \
        healpy \
        h5py \
        numexpr


USER $IMAGE_USER

COPY --chown=1000 *.sh $SRC/bin/
ENV SPT3G_SOFTWARE_PATH=/src/spt3g_software \
    SPT3G_SOFTWARE_BUILD_PATH=/src/spt3g_software/build
ENV PATH=${SPT3G_SOFTWARE_BUILD_PATH}/bin:$PATH \
    LD_LIBRARY_PATH=${SPT3G_SOFTWARE_BUILD_PATH}/spt3g:$LD_LIBRARY_PATH \
    PYTHONPATH=${SPT3G_SOFTWARE_BUILD_PATH}:$PYTHONPATH
