#!/bin/sh

[ -z "$OS_ARCH" ] && eval `/cvmfs/spt.opensciencegrid.org/py3-v1/setup.sh`
if [ -z "$SPT3G_SOFTWARE_PATH" ]; then
	export SPT3G_SOFTWARE_PATH=$(cd `dirname $0`/../../..; pwd)
	export PYTHONPATH=$SPT3G_SOFTWARE_PATH/../$OS_ARCH:$PYTHONPATH # XXX hard assumption about build layout!
fi

set -e

cd `dirname $0`/../../../chicago_processing

PYTHONPATH=.:$PYTHONPATH python scheduler.py /spt/user/nwhitehorn/sptpol/ /spt/user/nwhitehorn/sptpol/autoproc.txt /spt/user/nwhitehorn/sptpol/autoproc/ /spt/user/nwhitehorn/sptpol/ fullrate /tmp/x509up_u51232 $@

