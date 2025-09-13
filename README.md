# Cosmological Data Analysis Container

```bash
apptainer build --sandbox {container_path} docker://ghcr.io/cailmdaley/containers:latest
```

```bash
docker pull ghcr.io/cailmdaley/containers:latest
docker run --rm -it ghcr.io/cailmdaley/containers:latest bash
```

Tools included:
- Astronomy: `astropy`, `regions`, `reproject`, `skyproj`, `sip_tpv`, `sf_tools`, `lenspack`
- Observational cosmology: `healpy`, `healsparse`, `pymaster`, `pyccl`
  - CMB: `lenspyx`, `ducc0`, `camb`
  - Cosmic shear: `treecorr`, `pyccl`, `clmm`, `ngmix`, `shear_psf_leakage`
- Data analysis: `numpy`, `numba`, `scipy`, `pandas`, `h5py`, `pyarrow`, `joblib`, `seaborn`, `statsmodels`, `uncertainties`, `tqdm`
- Interactive: `ipython`, `jupyter`, `jupyterlab`, `jupytext`, `PyQt5`, `pyqtgraph`
- Writing tools: `pandoc`, `quarto`
- Low-level / compiler tools: `CMake`, `FFTW`, `CFITSIO`, `OpenMPI`, `gcc/g++/gfortran`
- Utilities: `vim`, `tmux`, `npm`, `xterm`, `rsync`
- Agent tools: `ripgrep`, `fd` (via `fd-find`), `fzf`, `bat`, `jq`, `yq`, `delta`, `httpie`, `aria2`, `sqlite3`, `hyperfine`, `btop`, `duf`, `procps`, `htop`, `less`, `unzip`, `zip`, `tldr`