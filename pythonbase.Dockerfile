FROM debian:buster-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential  \
        ca-certificates \
        cmake \
        curl \
        ffmpeg \
        gcc \
        git \
        gfortran \
        libbz2-dev \
        libffi-dev \
        liblzma-dev \
        libncurses5-dev \
        libreadline-dev \
        libssl-dev \
        libsqlite3-dev \
        make \
        man \
        nodejs \
        npm  \
        tar \
        tk-dev \
        vim \
        wget \
        xz-utils \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*


# setup an unprivileged user
ARG USER=cailmdaley
ENV USER $USER
ENV HOME /home/$USER
RUN adduser --disabled-password --gecos "Default user" $USER

# install python with pyenv since we need a dynamically-linked executable so
# that PyJulia works
ENV PATH="$HOME/.pyenv/shims:$HOME/.pyenv/bin:$PATH"
RUN curl https://pyenv.run | bash \
    && CFLAGS="-O2" PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.7.3 \
    && pyenv global 3.7.3

# install Python packages
RUN pip install --no-cache-dir \
        pip \
        jupyterlab \
        matplotlib \
        numpy \
        scipy \
    && jupyter labextension install @jupyterlab/toc \
    && rm -rf $HOME/.cache

USER $USER
WORKDIR $HOME
ENTRYPOINT ["/bin/bash"]
