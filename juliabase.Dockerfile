FROM debian:buster-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        cmake \
        curl \
        gcc \
        git \
        ffmpeg \
        gfortran \
        make \
        man \
        tar \
        vim \
        wget \
    && rm -rf /var/lib/apt/lists/*

# Install Julia
RUN mkdir /opt/julia \
    && curl -L https://julialang-s3.julialang.org/bin/linux/x64/1.3/julia-1.3.0-rc2-linux-x86_64.tar.gz | tar zxf - -C /opt/julia --strip=1 \
    && ln -s /opt/julia/bin/julia /usr/local/bin

# setup an unprivileged user
ARG USER=cailmdaley
ENV USER $USER
ENV HOME /home/$USER
RUN adduser --disabled-password --gecos "Default user" $USER
USER $USER
WORKDIR $HOME

# install Julia packages
ENV JULIA_PROJECT=$HOME
COPY --chown=1000 Project.toml $HOME/
RUN julia -e 'using Pkg; pkg"instantiate; precompile"'
RUN mkdir $HOME/src