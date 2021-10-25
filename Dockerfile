FROM python:3.9 as coscine

### Base Container
################################################################################

# install apt-get dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        curl \
        gcc \
        g++ \
        gfortran \
        git \
        less \
        libfontconfig \
        man \
        make \
        sudo \
        tar \
        vim \
        wget \
        xz-utils

# add very useful permissions tidier
COPY fix-permissions /usr/local/bin/fix-permissions
RUN chmod a+rx /usr/local/bin/fix-permissions

# create unprivileged user $IMAGE_USER
ENV IMAGE_USER=imageuser \
    IMAGE_UID=1000 \
    IMAGE_GID=100
ENV HOME=/home/$IMAGE_USER
RUN useradd -m -s /bin/bash -N -u $IMAGE_UID $IMAGE_USER && \
    fix-permissions $HOME
WORKDIR $HOME

# create user-writable /src and /src/bin directory for software and config files
ENV SRC=/src
RUN mkdir -p $SRC/bin \
    && chown -R $IMAGE_USER $SRC \
    && fix-permissions $SRC
ENV PATH=$SRC/bin:$PATH


### Python / Julia Setup
# all of the following can be done without sudo

# get julia binary
ENV JULIA_VERSION=1.7
COPY --from=julia:1.6 --chown=1000 /usr/local/julia /usr/local/julia
RUN ln -s /usr/local/julia/bin/julia /usr/local/bin/julia

USER $IMAGE_USER

# tell python where to --user install packages
ENV PYTHONUSERBASE=$SRC/.python
ENV PATH=$PYTHONUSERBASE/bin:$PATH

# install python and julia packages
RUN pip install --no-cache-dir \
        julia \
        matplotlib \
        numpy \
        pandas \
        scipy \
        snakemake

# remove default python command
CMD []


### Development Image: jupyter, writing tools
################################################################################
FROM coscine as coscine-dev

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        nodejs \
        npm \
        openssh-server

# Pandoc
ENV PANDOC_VERSION=2.14.1
RUN wget https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-amd64.deb \
    && dpkg -i pandoc-${PANDOC_VERSION}-1-amd64.deb \
    && rm pandoc-${PANDOC_VERSION}-1-amd64.deb
ENV XDG_DATA_HOME=$SRC

USER $IMAGE_USER

# TinyTeX
ENV TINYTEX_DIR=$SRC
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh"  | sh \
    && $SRC/.TinyTeX/bin/*/tlmgr option sys_bin $SRC/bin \
    && $SRC/.TinyTeX/bin/*/tlmgr path add \ 
    # install useful packages & pandoc pdf dependencies
    && tlmgr install \
        biblatex \
        caption \
        collection-fontsrecommended \
        dvipng \
        epsf \
        epstopdf-pkg \
        fancyhdr \
        grfext \
        lm-math \
        mathtools \
        microtype \
        physics \
        revtex4-1 \
        savesym \
        siunitx \
        textcase \
        ulem \
        unicode-math \
        wrapfig \
        xcolor

# Install jupyter
ENV JUPYTER_PATH=$SRC/.jupyter
ENV JUPYTER_CONFIG_DIR=$JUPYTER_PATH \
    JUPYTER_DATA_DIR=$JUPYTER_PATH
RUN pip install --no-cache-dir jupyterlab pytest

### SPTlab
################################################################################
FROM coscine-dev as sptlab

USER root
# dependencies for 3g-software & namaster
RUN apt-get update && apt-get install -y --no-install-recommends \
        libboost-all-dev \
        libfftw3-dev \
        libflac-dev \
        libgsl0-dev \
        libnetcdf-dev \
        libcfitsio-dev # for namaster

USER $IMAGE_USER
RUN pip install --no-cache-dir \
        astropy \
        camb \
        healpy \
        h5py \
        numexpr \
        pymaster
COPY --chown=1000 *.sh $SRC/bin/

# # Production container (precompiled Julia packages)
# ################################################################################
# FROM base as prod

# ENV JULIA_DEPOT_PATH=$SRC/.julia
# COPY --chown=1000 Project.toml $JULIA_DEPOT_PATH/environments/v${JULIA_VERSION}/Project.toml
# RUN git clone https://github.com/SouthPoleTelescope/CMBLensingSPT3GInterface /src/CMBLensingSPT3GInterface \
#     && julia -e 'using Pkg; Pkg.instantiate()' \
# #     && julia -e 'using Pkg; Pkg.develop(path="/src/CMBLensingSPT3GInterface")' \
#     && julia -e 'using Pkg; Pkg.precompile()'
