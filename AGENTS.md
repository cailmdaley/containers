# AGENTS.md

This file provides guidance to agents when working with code in this repository.

## Repository Overview

This repository contains a Docker-based cosmological data analysis environment designed for cosmic shear, CMB lensing, and cross-correlation analyses. The container provides a comprehensive scientific computing stack without conda, built on Python 3.11 with specialized astronomy and cosmology libraries.

## Development Workflow

### Making Changes
The workflow is straightforward:
1. Modify repository contents (typically `Dockerfile` or `scripts/`)
2. Commit changes
3. Push to trigger automated build

### Container Build Process
- **Build Location**: GitHub Actions only - no local builds or testing
- **Triggers**: Pushes to main/master/develop branches, tags, and manual dispatch
- **Pipeline**: Build → Smoke tests → Push to GitHub Container Registry (`ghcr.io/cailmdaley/containers`)
- **Testing**: Automated validation via `scripts/test.sh` during CI

## Architecture and Environment

### Container Design
- **Base**: Python 3.11 on Debian Bookworm
- **Build Strategy**: Single-stage image consolidating ShapePipe dependencies with agent-friendly CLI tools
- **Entry Point**: Bash shell with working directory at `/workspace`

### Scientific Computing Stack
The container is organized around three primary domains:

1. **Astronomy Libraries**: Core astronomical data handling (`astropy`, `regions`, `reproject`, `skyproj`)
2. **Observational Cosmology**: Specialized tools for cosmic shear analysis (`treecorr`, `pyccl`, `clmm`, `ngmix`) and CMB analysis (`lenspyx`, `ducc0`, `camb`)
3. **Data Analysis**: Standard scientific Python stack with performance libraries (`numpy`, `numba`, `scipy`, `pandas`)

### Developer Tools Integration
The container includes modern CLI tools optimized for agent workflows:
- **Search/Navigation**: `ripgrep`, `fd`, `fzf`, `bat`
- **Data Tools**: `jq`, `yq`, `sqlite3`, `httpie`
- **System Monitoring**: `btop`, `htop`, `hyperfine`
- **Text Processing**: `pandoc`, `quarto` for scientific writing

### Dependency Management
- **Git-based packages**: Several packages installed directly from Git repositories for latest features
- **No conda**: Pure pip-based installation for faster builds and simpler dependency resolution
- **Relaxed pins**: Dependencies use minimal version constraints for rapid development

## Key Files

- **`Dockerfile`**: Main container definition with all dependencies and tools
- **`scripts/test.sh`**: Smoke tests that validate core library imports and tool versions
- **`.github/workflows/deploy-image.yml`**: CI pipeline for automated building and publishing