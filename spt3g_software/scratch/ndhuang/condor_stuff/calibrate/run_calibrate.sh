#!/bin/sh

# Some directories
CODE=$_CONDOR_SCRATCH_DIR/code
mkdir $CODE
STARTDIR=$PWD
OUT=$_CONDOR_SCRATCH_DIR/output
mkdir $OUT

# Debug/status info
echo Host: `hostname`
echo OS/Arch: `uname -a`
echo OS_ARCH: $OS_ARCH

# set up environment
eval `/cvmfs/spt.opensciencegrid.org/py3-v1/setup.sh`
# Set up proxy if needed
# export connectHTTPProxy=${connectHTTPProxy:-UNAVAILABLE}
# if [ "$connectHTTPProxy" != UNAVAILABLE ]; then
#   export http_proxy=$connectHTTPProxy
# fi
# export OSG_SQUID_LOCATION=${OSG_SQUID_LOCATION:-UNAVAILABLE}
# if [ "$OSG_SQUID_LOCATION" != UNAVAILABLE ]; then
#   export http_proxy=$OSG_SQUID_LOCATION
# fi
echo OSG proxy:  $OSG_SQUID_LOCATION
echo Using proxy:  $http_proxy
export X509_USER_PROXY=$STARTDIR/x509up_u39371

set -e
# Get me some software
wget --quiet http://stash.ci-connect.net/spt/public/ndhuang/spt3g.tgz
# And set it up
tar -xf spt3g.tgz -C $CODE
export LD_LIBRARY_PATH=$CODE/spt3g/:$LD_LIBRARY_PATH
export PYTHONPATH=$CODE:$PYTHONPATH
# python -c 'from spt3g.util.genericutils import smooth, deriv, float_mode'


# Get me some datas
INPUT_URLS=""
LOCAL_PATH=""
for file in $1; do
    INPUT_URLS="$INPUT_URLS gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$file"
    LOCAL_INPUTS="$LOCAL_INPUTS data/`basename $file`"
done
mkdir data
for url in $INPUT_URLS; do
    echo Acquiring $url...
    globus-url-copy $url file://$_CONDOR_SCRATCH_DIR/data/
done

# Do things?
OUTFILE="$OUT/`basename $2`"
echo $OUTFILE
things="python $STARTDIR/analyze_calibrator.py $LOCAL_INPUTS -o $OUTFILE"
echo "Doing $things"
$things

# get datas back
echo "Sending $2 home"
globus-url-copy -cd file://$OUTFILE gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$2