#!/bin/sh

[ -z "$OS_ARCH" ] && eval `/software/clustertools/py3-v3/setup.sh`
if [ -z "$SPT3G_SOFTWARE_PATH" ]; then
	export SPT3G_SOFTWARE_PATH=$(cd `dirname $0`/..; pwd)
	export PYTHONPATH=$SPT3G_SOFTWARE_PATH/build_py3/$OS_ARCH:$PYTHONPATH
fi

set -e

cd `dirname $0`

LOGFILE=/poleanalysis/sptdaq/calresult/logs/scheduler.log
PYTHONPATH=.:$PYTHONPATH python -u scheduler.py /data/autoproc/db /data/autoproc/db/autoproc.txt /poleanalysis/sptdaq/calresult /spt_data/bolodata fullrate '' --verbose $@ >> $LOGFILE
