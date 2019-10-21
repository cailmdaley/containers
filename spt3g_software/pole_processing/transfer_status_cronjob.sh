#!/bin/sh

eval `/software/clustertools/py3-v3/setup.sh`
export PYTHONUNBUFFERED=1

ENVSH=$HOME/code/spt3g_software/build_py3/$OS_ARCH/env-shell.sh
PYDIR=$HOME/code/spt3g_software/pole_processing
$ENVSH python $PYDIR/check_sent.py
