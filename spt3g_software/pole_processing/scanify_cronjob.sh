#!/bin/sh

eval `/software/clustertools/py3-v3/setup.sh`
export PYTHONUNBUFFERED=1

ENVSH=$HOME/code/spt3g_software/build_py3/$OS_ARCH/env-shell.sh
PYDIR=$HOME/code/spt3g_software/pole_processing

sleep 20

# scanifier
LOGFILE=/data/autoproc/logs/scanify/scanify_submit.log
$ENVSH python $PYDIR/scanify_manager.py update --transfer-fullrate-cal --transfer-fullrate-fields 1 >> $LOGFILE
$ENVSH python $PYDIR/post_process_scanify.py --keep-after 12 >> $LOGFILE
