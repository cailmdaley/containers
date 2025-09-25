# syntax=docker/dockerfile:1.7

# Single-stage image: merges ShapePipe's Dockerfile with SP Validation additions
# and a rich CLI toolset for agents. No conda.

FROM python:3.11-bookworm

LABEL description="Cosmological data analysis environment: cosmic shear, CMB lensing, cross-correlations"

ENV SHELL=/bin/bash \
    DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    LANG=C.UTF-8

RUN set -eux; \
    # --- System upgrade --- \
    apt-get update -y --quiet --fix-missing && \
    apt-get dist-upgrade -y --quiet --fix-missing && \
    # --- Base build and runtime deps (ShapePipe core) --- \
    pkgs="\
      apt-utils \
      autoconf \
      automake \
      build-essential \
      chafa \
      cmake \
      curl \
      wget \
      ffmpeg \
      g++ \
      gcc \
      gfortran \
      git \
      git-lfs \
      libatlas-base-dev \
      libblas-dev \
      liblapack-dev \
      libcfitsio-dev \
      libhealpix-cxx-dev \
      libfftw3-bin \
      libfftw3-dev \
      libgl1-mesa-glx \
      libtool \
      libtool-bin \
      libtool-doc \
      libsharp-dev \
      locales \
      locate \
      make \
      openmpi-bin \
      libopenmpi-dev \
      pkg-config \
      protobuf-compiler \
      psfex \
      source-extractor \
      weightwatcher \
      vim \
      xterm \
      libgsl-dev \
      npm \
      tmux \
      ca-certificates \
      pandoc \
      # Agentic / developer tools 
      ripgrep \
      fd-find \
      fzf \
      bat \
      jq \
      yq \
      delta \
      httpie \
      aria2 \
      sqlite3 \
      hyperfine \
      btop \
      duf \
      procps \
      htop \
      less \
      unzip \
      zip \
      rsync \
      openssh-client \
      tldr "; \
    # --- Install --- \
    apt-get install -y --quiet --no-install-recommends $pkgs && \
    # convenience symlinks (Debian names) \
    if command -v batcat >/dev/null 2>&1; then ln -sf "$(command -v batcat)" /usr/local/bin/bat; fi && \
    if command -v fdfind  >/dev/null 2>&1; then ln -sf "$(command -v fdfind)"  /usr/local/bin/fd;  fi && \
    # cleanup \
    apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

# Install modern Neovim
RUN set -eux; \
    curl -fsSL -o /tmp/nvim.tar.gz https://github.com/neovim/neovim/releases/download/stable/nvim-linux-x86_64.tar.gz; \
    tar -C /usr/local -xzf /tmp/nvim.tar.gz; \
    ln -sf /usr/local/nvim-linux-x86_64/bin/nvim /usr/local/bin/nvim; \
    rm -f /tmp/nvim.tar.gz

# Chafa (static Linux x86_64 build)
RUN set -eux; \
    curl -fsSL -o /tmp/chafa.tar.gz https://hpjansson.org/chafa/releases/static/chafa-1.16.2-1-x86_64-linux-gnu.tar.gz; \
    tar -C /usr/local -xzf /tmp/chafa.tar.gz; \
    rm -f /tmp/chafa.tar.gz; \
    ln -sf /usr/local/chafa-1.16.2-1-x86_64-linux-gnu/chafa /usr/local/bin/chafa

# PolSpice (spherical power spectrum estimator)
RUN set -eux; \
    tmpdir="$(mktemp -d)"; \
    cd "$tmpdir"; \
    curl -fsSLO https://www2.iap.fr/users/hivon/software/PolSpice/ftp/PolSpice_v03-08-03.tar.gz; \
    tar -xzf PolSpice_v03-08-03.tar.gz; \
    cd PolSpice_v03-08-03; \
    mkdir build; \
    cd build; \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DHEALPIX=/usr; \
    cmake --build . --parallel "$(nproc)"; \
    install_dir=/opt/polspice; \
    mkdir -p "$install_dir"; \
    cp -r ../bin/. "$install_dir"/; \
    install -m 0755 "$install_dir/spice" /usr/local/bin/spice; \
    ln -sf spice /usr/local/bin/polspice; \
    rm -rf "$tmpdir"

# Zellij terminal multiplexer
COPY zellij-config.kdl /etc/zellij/config.kdl
RUN set -eux; \
    mkdir -p /tmp/zellij && \
    cd /tmp/zellij && \
    wget https://github.com/zellij-org/zellij/releases/download/v0.43.1/zellij-no-web-x86_64-unknown-linux-musl.tar.gz && \
    tar -xzf zellij-no-web-x86_64-unknown-linux-musl.tar.gz && \
    cp zellij /usr/local/bin/ && \
    rm -rf /tmp/zellij

# Quarto (for scientific writing) — prerelease channel, AMD64
RUN set -eux; \
    ARCH=amd64; \
    curl -L "$(curl -fsSL https://quarto.org/docs/download/_prerelease.json | grep -oP "(?<=\\\"download_url\\\":\\s\\\")https.*${ARCH}\\.deb")" -o /tmp/quarto.deb; \
    dpkg -i /tmp/quarto.deb; \
    rm /tmp/quarto.deb

# Python dependencies (relaxed pins for rapid development) + Snakemake
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir \
      astropy \
      adjustText \
      camb \
      clmm \
      colorama \
      cs_util \
      emcee \
      h5py \
      healpy \
      healsparse \
      importlib_metadata \
      ipython \
      joblib \
      jupyter \
      jupyterlab \
      jupytext \
      lenspack \
      lmfit \
      numexpr \
      numpy \
      numba \
      mpi4py \
      numpydoc \
      opencv-python-headless \
      pandas \
      PyQt5 \
      pyqtgraph \
      pyccl \
      pyarrow \
      pymaster \
      pytest \
      pytest-cov \
      pytest-pydocstyle \
      pytest-pycodestyle \
      regions \
      reproject \
      scipy \
      seaborn \
      skyproj \
      sf_tools \
      sip_tpv \
      skaha \
      sqlitedict \
      statsmodels \
      termcolor \
      tqdm \
      treecorr \
      uncertainties \
      vos \
      snakemake \
      snakemake-executor-plugin-slurm \
      lenspyx \
      ducc0 \
      # Git-based dependencies
      git+https://github.com/aguinot/cosmo-numba \
      git+https://github.com/benabed/getdist.git@upper_triangle_whisker \
      git+https://github.com/CosmoStat/shear_psf_leakage.git@develop \
      git+https://github.com/aguinot/ngmix@stable_version \
      git+https://github.com/tobias-liaudat/Stile@v0.1 \
      git+https://github.com/CosmoStat/sp_validation

# Working directory for user code
WORKDIR /workspace

# Default entry
ENTRYPOINT ["bash"]
CMD ["-l"]
