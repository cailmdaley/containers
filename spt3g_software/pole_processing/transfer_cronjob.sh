#!/bin/sh

eval `/software/clustertools/py3-v3/setup.sh`
export PYTHONUNBUFFERED=1

if [ -n "$SPT3G_SOFTWARE_PATH" ]; then
    unset SPT3G_SOFTWARE_PATH
fi

ENVSH=$HOME/code/spt3g_software/build_py3/$OS_ARCH/env-shell.sh
PYDIR=$HOME/code/spt3g_software/pole_processing
LOGDIR=/data/autoproc/logs/transfer

sleep 30

# observation transfer
LOGFILE=$LOGDIR/transfer.log
$ENVSH python $PYDIR/transfer_manager.py >> $LOGFILE &

sleep 10

# aux file transfer
LOGFILE=$LOGDIR/aux_transfer.log
$ENVSH python $PYDIR/aux_transfer_manager.py >> $LOGFILE &

sleep 10

# verification
LOGFILE=$LOGDIR/verify.log
$ENVSH python $PYDIR/verify_transfer.py --logdir $LOGDIR >> $LOGFILE &

sleep 10

# rsync
LOGFILE=$LOGDIR/rsync.log
$ENVSH python $PYDIR/rsync_manager.py >> $LOGFILE &

wait
