#!/bin/bash

# NB: this script must be run on {scott,amundsen}.grid.uchicago.edu

eval `/cvmfs/spt.opensciencegrid.org/py3-v3/setup.sh`
export PYTHONUNBUFFERED=1

PYDIR=$HOME/code/spt3g_software_stable/pole_processing

$PYDIR/unpack_pydfmux 2>&1 | mail -E -s "Cron $HOSTNAME: Pydfmux Archive Extraction Summary" $LOGNAME
$PYDIR/unpack_tar 2>&1 | mail -E -s "Cron $HOSTNAME: User Archive Extraction Summary" $LOGNAME
