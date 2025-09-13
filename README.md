# Cosmological Data Analysis Container

Apptainer sandbox build (copy/paste):

```bash
apptainer build --sandbox /n17data/cdaley/containers/containers docker://ghcr.io/cailmdaley/containers:latest
```

Key tool themes:
- Astronomy: `astropy`, `regions`, `reproject`, `skyproj`, `sip_tpv`, `sf_tools`, `lenspack`
- Observational cosmology: `healpy`,  `healsparse`, `pymaster`, `pyccl`
  - CMB: `lenspyx`, `ducc0`, `camb`
  - Cosmic shear: `treecorr`, `pyccl`, `clmm`, `ngmix`, `shear_psf_leakage`
- Data analysis: `numpy`, `numba`, `scipy`, `pandas`, `h5py`, `pyarrow`, `joblib`, `seaborn`, `statsmodels`, `uncertainties`, `tqdm`
- Interactive: `ipython`, `jupyter`, `jupyterlab`, `jupytext`, `PyQt5`, `pyqtgraph`
- Writing tools: `pandoc`, `quarto`
- Low-level / compiler tools: `CMake`, `FFTW`, `CFITSIO`, `OpenMPI`, `gcc/g++/gfortran`
- Utilities: `vim`, `tmux`, `npm`, `xterm`, `rsync`
- Agent tools: `ripgrep`, `fd` (via `fd-find`), `fzf`, `bat`, `jq`, `yq`, `delta`, `httpie`, `aria2`, `sqlite3`, `hyperfine`, `btop`, `duf`, `procps`, `htop`, `less`, `unzip`, `zip`, `tldr`

## Build & Publish

This image is built exclusively by GitHub Actions and pushed to GHCR on every push to the main branches and on tags (see `.github/workflows/deploy-image.yml`). Use the workflow’s "Run workflow" button for manual builds.

## Use

Pull and run from GHCR:

```bash
docker pull ghcr.io/cailmdaley/containers:latest
docker run --rm -it ghcr.io/cailmdaley/containers:latest bash
```

## What’s Inside
- Base: `python:3.11-bookworm`.
- ShapePipe OS deps: FFTW, CFITSIO, OpenMPI, PSFEx, Source Extractor, WeightWatcher, toolchain, etc.
- Python stack (unpinned for rapid dev): full ShapePipe stack plus all `sp_validation` dependencies, e.g., astropy, galsim, treecorr, numpy/numba, healpy/healsparse, pyccl, pymaster, camb, clmm, jupyter/jupyterlab, pytest, etc. Also includes Lenspyx and ducc0 for CMB lensing.
- Snakemake (from SP Validation).
- Agentic/dev tools (Debian): ripgrep, fd-find→`fd`, fzf, bat→`bat`, jq, yq, delta, httpie, aria2, sqlite3, hyperfine, btop, duf, procps, htop, less, unzip/zip, rsync, openssh-client, tldr.

All tools are included by default to empower automated agents and interactive workflows.
