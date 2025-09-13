#!/usr/bin/env bash
set -euo pipefail
echo "[smoke] Python imports"
python - <<'PY'
import importlib, sys
mods = [
    'astropy', 'healpy', 'healsparse', 'pyccl', 'pymaster', 'camb', 'clmm',
    'treecorr', 'lenspyx', 'ducc0', 'numpy', 'numba', 'mpi4py', 'sp_validation'
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

echo "[smoke] pandoc --version"
pandoc --version | head -n 1

echo "[smoke] quarto --version"
quarto --version

echo "Smoke tests passed"
