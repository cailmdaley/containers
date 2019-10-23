FROM cailmdaley/juliabase:latest

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential  \
        libbz2-dev       \
        libffi-dev       \
        liblzma-dev      \
        libncurses5-dev  \
        libreadline-dev  \
        libssl-dev       \
        libsqlite3-dev   \
        nodejs           \
        npm              \
        tk-dev           \
        xz-utils         \
        zlib1g-dev       \
    && rm -rf /var/lib/apt/lists/*

# install python with pyenv since we need a dynamically-linked executable so
# that PyJulia works
ENV PATH="$HOME/.pyenv/shims:$HOME/.pyenv/bin:$PATH"
user $IMAGEUSER
RUN curl https://pyenv.run | bash \
    && CFLAGS="-O2" PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.7.3 \
    && pyenv global 3.7.3 \
    && chmod -R 777 $HOME

# install Python packages
RUN pip install --no-cache-dir \
        pip \
        jupyterlab \
        matplotlib \
        numpy \
        scipy \
        git+https://github.com/JuliaPy/pyjulia.git \
    && jupyter labextension install @jupyterlab/toc \
    && rm -rf $HOME/.cache

# install Julia packages
RUN julia -e 'using Pkg; pkg"add PyCall PyPlot; precompile"' \
    && chmod -R 777 $HOME
