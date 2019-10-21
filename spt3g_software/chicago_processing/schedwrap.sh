#!/bin/sh

[ -z "$OS_ARCH" ] && eval `/cvmfs/spt.opensciencegrid.org/py3-v2/setup.sh`
if [ -z "$SPT3G_SOFTWARE_PATH" ]; then
	export SPT3G_SOFTWARE_PATH=$(cd `dirname $0`/..; pwd)
	export PYTHONPATH=$SPT3G_SOFTWARE_PATH/../$OS_ARCH:$PYTHONPATH # XXX hard assumption about build layout!
fi

set -e

cd `dirname $0`

# Issue proxy cert for file transfer
[ -d proxycerts ] || mkdir proxycerts
X509_USER_CERT=~/.globus/autoproc_unl.pem X509_USER_KEY=~/.globus/autoproc_unl_key.pem grid-proxy-init -q -pwstdin -out proxycerts/proxycert_$$.p12 -verify < /dev/null > /dev/null

PYTHONPATH=.:$PYTHONPATH python scheduler.py --verbose /spt/data/transfer_database/ /spt/user/production/autoproc.txt /spt/user/production/ /spt/data/bolodata/ fullrate proxycerts/proxycert_$$.p12 $@
