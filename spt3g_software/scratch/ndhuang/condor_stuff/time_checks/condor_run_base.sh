#!/bin/sh
# $1 input
# $2 output
# $3 args to script

# Some directories
CODE=$_CONDOR_SCRATCH_DIR/code
STARTDIR=$PWD

# Debug/status info
echo Host: `hostname`
echo OS/Arch: `uname -a`
echo OS_ARCH: $OS_ARCH

# set up environment
eval `/cvmfs/spt.opensciencegrid.org/py3-v1/setup.sh`
# Set up proxy if needed
export connectHTTPProxy=${connectHTTPProxy:-UNAVAILABLE}
if [ "$connectHTTPProxy" != UNAVAILABLE ]; then
  export http_proxy=$connectHTTPProxy
fi
export OSG_SQUID_LOCATION=${OSG_SQUID_LOCATION:-UNAVAILABLE}
if [ "$OSG_SQUID_LOCATION" != UNAVAILABLE ]; then
  export http_proxy=$OSG_SQUID_LOCATION
fi
echo OSG proxy:  $OSG_SQUID_LOCATION
echo Using proxy:  $http_proxy
export X509_USER_PROXY=$STARTDIR/x509up_u39371

# Get me some software
wget http://stash.osg-connect.net/spt/public/ndhuang/spt3g.tgz
# And set it up
tar xf spt3g.tgz -C $CODE
$CODE/env-shell.sh

# Get me some datas
INPUT_URLS=""
LOCAL_PATH=""
for file in $2; do
        INPUT_URLS="$INPUT_URLS gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$file"
        LOCAL_PATH="$LOCAL_PATH data/`basename $file`"
done
mkdir data
for url in $INPUT_URLS; do
        echo Acquiring $url...
        globus-url-copy $url file://$_CONDOR_SCRATCH_DIR/data/
done

# Do things?

# get datas back
globus-url-copy -cd `basename $4` gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$4