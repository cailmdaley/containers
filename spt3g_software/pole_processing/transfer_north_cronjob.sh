#!/bin/sh

# NB: this script must be run as user spt on spt-buffer.grid.uchicago.edu

eval `/cvmfs/spt.opensciencegrid.org/py3-v3/setup.sh`
export PYTHONUNBUFFERED=1

ENVSH=$HOME/git/spt3g_software/build/env-shell.sh
PYDIR=$HOME/git/spt3g_software/pole_processing
LOGDIR=/buffer/transfer_logs

# archive
LOGFILE=$LOGDIR/archive_north.log
$ENVSH python $PYDIR/archive_north.py >> $LOGFILE &

sleep 5

# observation copy
LOGFILE=$LOGDIR/transfer_north_copy.log
$ENVSH python $PYDIR/transfer_manager_north_copy.py >> $LOGFILE &

sleep 5

# aux file copy
LOGFILE=$LOGDIR/aux_transfer_north_copy.log
$ENVSH python $PYDIR/aux_transfer_manager_north_copy.py >> $LOGFILE &

sleep 5

# observation transfer
LOGFILE=$LOGDIR/transfer_north.log
$ENVSH python $PYDIR/transfer_manager_north.py >> $LOGFILE &

sleep 5

# rsync file transfer
LOGFILE=$LOGDIR/rsync_north.log
$ENVSH python $PYDIR/rsync_manager_north.py >> $LOGFILE &

sleep 5

# verification
LOGFILE=$LOGDIR/verify_north.log
$ENVSH python $PYDIR/verify_transfer_north.py >> $LOGFILE &

wait

OLDFILES=`find /buffer/3g_verified/ -maxdepth 1 -type f -regextype posix-extended -regex '/buffer/3g_verified/[0-9]{8}\.tar' -mtime +15 -print`
if [ -n "$OLDFILES" ]; then
    echo "Removing old data tarfiles from buffer"
    rm -v $OLDFILES
fi
