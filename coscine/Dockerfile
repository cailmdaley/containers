FROM python:3.9

### Basic Setup
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
        ssh-client \
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
################################################################################

# all of the following can be done without sudo

# get julia binary
# ENV JULIA_DEPOT_PATH=$SRC/.julia
COPY --from=julia:1.6.0 --chown=1000 /usr/local/julia /usr/local/julia
RUN ln -s /usr/local/julia/bin/julia /usr/local/bin/julia

USER $IMAGE_USER

# tell python where to --user install packages
ENV PYTHONUSERBASE=$SRC/.python
ENV PATH=$PYTHONUSERBASE/bin:$PATH

# tell jupyter where to put things
ENV JUPYTER_PATH=$SRC/.jupyter
ENV JUPYTER_CONFIG_DIR=$JUPYTER_PATH \
    JUPYTER_DATA_DIR=$JUPYTER_PATH

# install python and julia packages
RUN pip install --no-cache-dir \
        julia \
        jupyterlab \
        matplotlib \
        numpy \
        scipy
    # && julia -e 'using Pkg; pkg"add PyCall IJulia; precompile"'
