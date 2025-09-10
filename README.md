# Cosmological Data Analysis Container

Single Dockerfile, no conda. Merges the ShapePipe environment with SP Validation dependencies, adds agent-friendly CLI tools, and targets cosmic shear, CMB lensing, and cross-correlations. No multi-stage builds.

## Build & Publish

This image is built exclusively by GitHub Actions and pushed to GHCR on every push to the main branches and on tags (see `.github/workflows/deploy-image.yml`). Use the workflow’s "Run workflow" button for manual builds.

## Use

Pull and run from GHCR:

```bash
docker pull ghcr.io/cailmdaley/containers:latest
docker run --rm -it ghcr.io/cailmdaley/containers:latest bash
```

## What’s Inside
- Base: `python:3.12-bookworm`.
- ShapePipe OS deps: FFTW, CFITSIO, OpenMPI, PSFEx, Source Extractor, WeightWatcher, toolchain, etc.
- Python stack (unpinned for rapid dev): full ShapePipe stack plus all `sp_validation` dependencies, e.g., astropy, galsim, treecorr, numpy/numba, healpy/healsparse, pyccl, pymaster, camb, clmm, jupyter/jupyterlab, pytest, etc. Also includes Lenspyx and ducc0 for CMB lensing.
- Snakemake (from SP Validation).
- Agentic/dev tools (Debian): ripgrep, fd-find→`fd`, fzf, bat→`bat`, jq, yq, delta, httpie, aria2, sqlite3, hyperfine, btop, duf, procps, htop, less, unzip/zip, rsync, openssh-client, tldr.

All tools are included by default to empower automated agents and interactive workflows.

## HPC via Apptainer (optional)

```bash
# pull from GHCR (recommended)
apptainer pull cosmology-env.sif docker://ghcr.io/cailmdaley/containers:latest
apptainer exec cosmology-env.sif bash
```

Or publish to a registry (e.g. GHCR) and pull with `apptainer pull docker://<registry>/<image>:<tag>`.
