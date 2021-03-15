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


# Full development container
################################################################################
FROM base as dev
RUN pip install --no-cache-dir pytest


# Production container (precompiled Julia packages)
################################################################################
FROM base as prod

ENV JULIA_DEPOT_PATH=$SRC/.julia
COPY --chown=1000 Project.toml $JULIA_DEPOT_PATH/environments/v${JULIA_VERSION}/Project.toml
RUN git clone https://github.com/SouthPoleTelescope/CMBLensingSPT3GInterface /src/CMBLensingSPT3GInterface \
    && julia -e 'using Pkg; Pkg.instantiate()' \
#     && julia -e 'using Pkg; Pkg.develop(path="/src/CMBLensingSPT3GInterface")' \
    && julia -e 'using Pkg; Pkg.precompile()'


