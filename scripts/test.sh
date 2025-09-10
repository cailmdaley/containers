#!/usr/bin/env bash
set -euo pipefail
echo "[smoke] Python imports"
python - <<'PY'
import importlib, sys
mods = [
    'astropy', 'healpy', 'healsparse', 'pyccl', 'pymaster', 'camb', 'clmm',
    'treecorr', 'lenspyx', 'ducc0', 'numpy', 'numba', 'mpi4py'
]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as e:
        bad.append((m, str(e)))
if bad:
    for m, err in bad:
        print(f"FAIL import {m}: {err}")
    sys.exit(1)
print("OK: core libraries import")
PY

echo "[smoke] snakemake --version"
snakemake --version

echo "[smoke] treecorr --help"
treecorr --help >/dev/null

echo "[smoke] healpy version"
python - <<'PY'
import healpy as hp
print("healpy", hp.__version__)
PY

echo "Smoke tests passed"
