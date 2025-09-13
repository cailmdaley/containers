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
      libfftw3-bin \
      libfftw3-dev \
      libgl1-mesa-glx \
      libtool \
      libtool-bin \
      libtool-doc \
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

# Quarto (for scientific writing) — latest AMD64 build
RUN set -eux; \
    STABLE_URL="https://github.com/quarto-dev/quarto-cli/releases/latest/download/quarto-linux-amd64.deb"; \
    TMP_DEB="/tmp/quarto-linux-amd64.deb"; \
    if ! curl -fL --retry 5 --retry-all-errors -o "$TMP_DEB" "$STABLE_URL"; then \
      echo "Stable Quarto download failed; trying prerelease channel"; \
      PRE_URL=$(curl -fsSL https://quarto.org/docs/download/_prerelease.json \
        | sed -n 's/.*"download_url"[[:space:]]*:[[:space:]]*"\([^"']*amd64\.deb\)".*/\1/p' \
        | head -n1); \
      test -n "$PRE_URL"; \
      curl -fL --retry 5 --retry-all-errors -o "$TMP_DEB" "$PRE_URL"; \
    fi; \
    apt-get update -y --quiet; \
    apt-get install -y --quiet --no-install-recommends "$TMP_DEB"; \
    rm -f "$TMP_DEB" && \
    apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

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
