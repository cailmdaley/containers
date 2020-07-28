FROM cailmdaley/coscine-tweaked:latest

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        git-lfs \
        g++ \
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


COPY --chown=1000 github-key $HOME/.ssh/id_rsa
RUN ssh-keyscan -t rsa github.com >> ~/.ssh/known_hosts \
    && git -C $SRC clone git@github.com:SouthPoleTelescope/spt3g_software.git \
    && cd $SRC/spt3g_software \
        && mkdir build \
        && cd build \
        && cmake .. \
        && make -j8

USER $IMAGE_USER
RUN $SRC/spt3g_software/build/env-shell.sh
ENV SPT3G_SOFTWARE_PATH=/src/spt3g_software \
	SPT3G_SOFTWARE_BUILD_PATH=/src/spt3g_software/build
ENV PATH=${SPT3G_SOFTWARE_BUILD_PATH}/bin:$PATH \
    LD_LIBRARY_PATH=${SPT3G_SOFTWARE_BUILD_PATH}/spt3g:$LD_LIBRARY_PATH \
    PYTHONPATH=${SPT3G_SOFTWARE_BUILD_PATH}:$PYTHONPATH
