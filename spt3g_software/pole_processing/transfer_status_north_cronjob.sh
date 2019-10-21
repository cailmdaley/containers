#!/bin/sh

# NB: this script must be run as user spt on spt-buffer.grid.uchicago.edu

eval `/cvmfs/spt.opensciencegrid.org/py3-v3/setup.sh`
export PYTHONUNBUFFERED=1

ENVSH=$HOME/git/spt3g_software/build/env-shell.sh
PYDIR=$HOME/git/spt3g_software/pole_processing
$ENVSH python $PYDIR/check_sent.py | mail -E -s "SPTR Transfer Status" -c tcrawfor@kicp.uchicago.edu $LOGNAME
